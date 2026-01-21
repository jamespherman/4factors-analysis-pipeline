# DATA_SPEC.md - Data Format and Export Specification

## Overview

This document specifies the data format for RNN fitting and the MATLAB export script that generates it from the existing 4factors-analysis-pipeline.

## Source Data

**Location**: `OneDrive/Neuronal Data Analysis/4factors-analysis-pipeline/data/processed/`

**Initial Session**: `Feynman_08_15_2025_SC/session_data.mat`

The MATLAB pipeline produces `session_data.mat` files containing:
- `trialInfo`: Trial metadata (target location, reward, salience, identity, etc.)
- `eventTimes`: Timestamps for task events (fixOn, targetOn, fixOff, saccadeOnset, reward)
- `spikes`: Spike times for each neuron
- `analysis`: Preprocessing outputs including neuron classification

## Exported Data Format

The MATLAB export script produces a single `.mat` file per session with the following structure:

### `rnn_export_<session_name>.mat`

```matlab
% Neural data
firing_rates       % [nNeurons × nTimeBins × nTrials] - Binned firing rates (sp/s)
neuron_ids         % [nNeurons × 1] - Original neuron indices from session_data
neuron_type        % [nNeurons × 1] - 1 = classic/excitatory, 2 = putative_interneuron/inhibitory
brain_area         % [nNeurons × 1] - Brain area code (1 = SC for this project)

% Trial information
n_trials           % scalar - Number of trials
n_neurons          % scalar - Number of neurons
n_time_bins        % scalar - Number of time bins per trial

% Time parameters
bin_size_ms        % scalar - Bin size in milliseconds (25)
time_axis          % [nTimeBins × 1] - Time relative to fixation onset (ms)
                   %   Note: bin 0 starts at (fixation_onset - pre_fix_ms)
                   %   So fixation occurs at bin index = pre_fix_ms / bin_size_ms
trial_duration_ms  % [nTrials × 1] - Actual duration of each trial (optional)

% Task inputs (time-varying) - ALL REQUIRED
input_fixation_on  % [nTimeBins × nTrials] - Binary, 1 when fixation point visible
input_target_loc   % [4 × nTimeBins × nTrials] - One-hot encoding of target location (1-4)
input_go_signal    % [nTimeBins × nTrials] - Binary, 1 after fixation offset
input_reward_on    % [nTimeBins × nTrials] - Binary, 1 during juice delivery

% Eye position (OPTIONAL - set to zeros or NaN if unavailable)
input_eye_x        % [nTimeBins × nTrials] - Eye position X (degrees), z-score normalized in Python
input_eye_y        % [nTimeBins × nTrials] - Eye position Y (degrees), z-score normalized in Python
                   % Note: If eye traces aren't available, fill with zeros. Model learns to ignore.

% Target features (revealed at target onset) - ALL REQUIRED
input_is_face      % [nTimeBins × nTrials] - Binary, 1 if face target (when target visible)
input_is_nonface   % [nTimeBins × nTrials] - Binary, 1 if non-face target
input_is_bullseye  % [nTimeBins × nTrials] - Binary, 1 if bullseye target
input_high_salience% [nTimeBins × nTrials] - Binary, 1 if high salience (bullseye only, 0 for face/nonface)
input_low_salience % [nTimeBins × nTrials] - Binary, 1 if low salience (bullseye only, 0 for face/nonface)

% Trial-level labels (for analysis, not model input)
trial_reward       % [nTrials × 1] - 1 = high, 0 = low
trial_probability  % [nTrials × 1] - 1 = high, 0 = low (optional)
trial_identity     % [nTrials × 1] - 1 = face, 2 = nonface, 3 = bullseye
                   %   Matches MATLAB stimType where face=1, nonface=2, bullseye>2
trial_salience     % [nTrials × 1] - 0 = N/A (face/nonface trials), 1 = high, 2 = low
                   %   Only bullseye trials have salience; face/nonface get 0
trial_location     % [nTrials × 1] - Target location index (1-4)
                   %   Matches targetLocIdx from MATLAB (1=loc1, 2=loc2, etc.)

% Metadata
session_name       % string - Session identifier
export_date        % string - Date of export
pipeline_version   % string - Version of analysis pipeline
```

## Temporal Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Time window start | -200ms pre-fixation | Capture baseline |
| Time window end | +500ms post-reward | Capture reward response |
| Bin size | 25ms | Matches existing pipeline; ~40 bins/sec |
| Alignment | Trial start (fixation onset) | All trials aligned to same reference |

## Trial Selection

**Include ALL successfully completed gSac_4factors trials** regardless of:
- Target laterality (all 4 locations included)
- Reward level
- Target identity/salience

**Rationale**: Putative interneurons may have bilateral or unusual spatial tuning. Location becomes an input to the network rather than a selection criterion.

## Neuron Classification Criteria

### Classic SC Neurons → Assigned to Excitatory Units

ALL criteria must be met:
1. **Task-responsive**: Significant increase (p<0.05) in visual OR delay OR saccade epoch vs baseline for ≥1 contralateral location
2. **Adequate firing**: Mean FR > 5 sp/s across session
3. **Unilateral**: NOT significantly activated at opposing locations (e.g., not both 135° and 315°)
4. **Not sparse**: ≤70% of 100ms bins empty across session

### Putative Interneurons → Assigned to Inhibitory Units

ALL criteria must be met:
1. **Task-modulated**: Significant effect somewhere in task (any epoch, any factor)
2. **Adequate firing**: Mean FR > 5 sp/s
3. **Fails classic**: Does not meet classic SC neuron criteria (suppression, bilateral, idiosyncratic)

### Neurons Excluded from RNN

- Mean FR < 5 sp/s (insufficient data)
- No task modulation (not informative)
- >70% empty bins (too sparse for reliable estimation)

## MATLAB Export Script

**Location**: `4factors-analysis-pipeline/scripts/export_for_rnn.m`

**Usage**:
```matlab
% Export single session
export_for_rnn('Feynman_08_15_2025_SC', output_path);

% Export with custom parameters
opts.bin_size_ms = 25;
opts.pre_fix_ms = 200;
opts.post_reward_ms = 500;
export_for_rnn(session_name, output_path, opts);
```

**Implementation Notes**:

1. Load `session_data.mat` from processed data directory
2. Apply neuron classification using existing `analysis.neuron_screening` results
3. Bin spike times into firing rates using specified bin size
4. Construct input signals from `eventTimes` and `trialInfo`
5. Handle variable trial lengths by:
   - Computing max trial duration across all trials
   - Padding shorter trials with NaN
   - Storing actual duration per trial for masking
6. Save all fields to single `.mat` file

## Python Data Loading

```python
import scipy.io as sio
import numpy as np

def load_rnn_data(filepath):
    """Load exported .mat file for RNN fitting.
    
    Returns:
        dict with numpy arrays for all fields
    """
    data = sio.loadmat(filepath, squeeze_me=True)
    
    # Remove MATLAB metadata fields
    data = {k: v for k, v in data.items() if not k.startswith('__')}
    
    # Verify shapes
    n_neurons = data['n_neurons']
    n_bins = data['n_time_bins']
    n_trials = data['n_trials']
    
    assert data['firing_rates'].shape == (n_neurons, n_bins, n_trials)
    
    return data
```

## Validation Checks

Before fitting, verify:

1. **No NaN in inputs**: All input signals should be valid (NaN only in padded regions)
2. **Firing rate sanity**: Mean FR should match expected range (1-50 sp/s typical for SC)
3. **E/I counts**: Should have both excitatory and inhibitory neurons classified
4. **Trial balance**: Check distribution across conditions (reward × location × identity)
5. **Temporal alignment**: Verify event times fall within expected windows

```python
def validate_export(data):
    """Run sanity checks on exported data."""
    
    # Check for NaN in critical fields
    assert not np.any(np.isnan(data['input_fixation_on']))
    
    # Check firing rate range
    mean_fr = np.nanmean(data['firing_rates'])
    assert 1 < mean_fr < 50, f"Suspicious mean FR: {mean_fr}"
    
    # Check neuron types
    n_exc = np.sum(data['neuron_type'] == 1)
    n_inh = np.sum(data['neuron_type'] == 2)
    print(f"Excitatory: {n_exc}, Inhibitory: {n_inh}")
    assert n_exc > 0 and n_inh > 0, "Need both E and I neurons"
    
    # Check trial counts per condition
    print(f"High reward trials: {np.sum(data['trial_reward'] == 1)}")
    print(f"Low reward trials: {np.sum(data['trial_reward'] == 0)}")
    
    return True
```

## File Naming Convention

```
rnn_export_<monkey>_<date>_<area>.mat

Examples:
- rnn_export_Feynman_08_15_2025_SC.mat
- rnn_export_Vito_03_22_2024_SC.mat
```

## Data Location

Exported files should be placed in:
```
4factors-rnn-analysis/data/
```

This directory is NOT tracked in git (add to .gitignore). The README in that directory should specify the path to the source data on the local machine.
