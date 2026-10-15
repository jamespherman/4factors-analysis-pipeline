# Neuron Screening Improvement Plan

## Overview

This document outlines a plan to improve the neuron screening pipeline with two main goals:
1. **Enable rapid threshold experimentation** by caching expensive computations (Friedman tests, task modulation)
2. **Improve screening accuracy** by detecting stable recording periods for each neuron and restricting analyses to those periods

## Background and Motivation

### Current Limitations

The current screening pipeline has two main limitations:

1. **Slow threshold iteration**: When experimenting with firing rate (FR) and sparsity thresholds, the pipeline recomputes Friedman tests (~23 min for SC sessions) even though these results are independent of the threshold values.

2. **Biased FR/sparsity estimates**: Some neurons exhibit idiosyncratic recording patterns:
   - A handful of spikes early in the session, then silence, then consistent firing
   - Gradual drift causing the neuron to "appear" or "disappear" mid-recording
   - Abrupt changes due to electrode movement or tissue settling

   Using first-spike-to-last-spike as the "active period" dramatically underestimates FR for neurons that were only stable for a subset of trials.

### Literature Review

Research on spike sorting quality metrics and recording stability provides several relevant approaches:

**Presence Ratio** ([SpikeInterface documentation](https://spikeinterface.readthedocs.io/en/latest/modules/qualitymetrics.html)):
- Defined as the fraction of time bins in which a neuron fires above a threshold
- Expected value for complete units: ≥90%
- Implementation: Divide recording into bins (e.g., 60s), count bins with ≥1 spike
- Low presence ratio indicates incompleteness or highly selective firing

**Amplitude Stability** ([Hill et al., 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3123734/)):
- Plot spike amplitude vs. time to detect drift
- Systematic changes indicate electrode movement or neuron state changes
- Flag units with gradual drifts or sudden jumps

**Drift Correction** ([MEDiCINe, eNeuro 2025](https://www.eneuro.org/content/12/3/ENEURO.0529-24.2025)):
- Brain tissue moves relative to electrodes due to pulsation, pressure changes, licking
- Modern approaches (Kilosort4, DREDge) estimate and correct for motion
- Neurons can be "lost" when they drift out of detection range

**Sliding Window Approaches** ([SpikeSift, 2025](https://arxiv.org/html/2504.01604v1)):
- Segment recordings into shorter intervals with minimal drift-induced variation
- Dynamically determine segment boundaries by analyzing amplitude fluctuations
- More robust than fixed user-defined windows

---

## Phase 1: Caching Friedman Test Results

### Objective
Cache Friedman test results to enable rapid iteration on FR/sparsity thresholds without recomputing statistical tests.

### Rationale
- Friedman tests depend only on spike timing relative to task epochs
- Results are independent of FR threshold and sparsity threshold
- Current runtime: ~200s per SC session × 7 sessions = ~23 min
- After caching: ~0s (load from disk)

### Implementation

#### 1.1 Define Cache Structure

Add to `session_data`:
```matlab
session_data.cached_friedman_results = struct(...
    'cluster_ids',        [],       % [nClusters × 1] for validation
    'sig_by_location',    [],       % [nClusters × nGroups × 3] boolean
    'group_mean_frs',     {},       % {1 × nGroups} cell array
    'location_angles',    [],       % [nGroups × 1] angles
    'epoch_frs',          [],       % [nClusters × 4 × nTrials] raw data
    'computation_params', struct(), % Parameters used (epoch windows, etc.)
    'computation_date',   ''        % Timestamp for reference
);
```

#### 1.2 Modify `screen_sc_neurons.m`

**Before Friedman test loop (around line 400):**
```matlab
% Check for cached Friedman results
use_cached_friedman = false;
if isfield(session_data, 'cached_friedman_results') && ...
        ~isempty(session_data.cached_friedman_results)
    cached = session_data.cached_friedman_results;
    % Validate: same clusters, same order
    if isequal(cached.cluster_ids(:), cluster_ids(:))
        use_cached_friedman = true;
        fprintf('screen_sc_neurons: Using cached Friedman test results.\n');
        sig_by_location = cached.sig_by_location;
        group_mean_frs = cached.group_mean_frs;
    end
end

if ~use_cached_friedman
    % Existing Friedman test computation (lines 412-452)
    % ... existing code ...

    % Cache results in session_data
    session_data.cached_friedman_results = struct(...
        'cluster_ids', cluster_ids, ...
        'sig_by_location', sig_by_location, ...
        'group_mean_frs', {group_mean_frs}, ...
        'location_angles', location_angles, ...
        'computation_date', datestr(now));
end
```

#### 1.3 Update Function Signature

Modify `screen_sc_neurons.m` to return `session_data`:
```matlab
function [selected_neurons, screening_info, neuron_class, session_data] = ...
    screen_sc_neurons(session_data, core_data, condition_defs, project_root)
```

#### 1.4 Update Caller (`apply_neuron_screening.m`)

Update call at line ~313 to capture returned `session_data`.

### Testing
1. Run screening with `force_rerun.screening = true`
2. Verify "Using cached Friedman test results" appears on second run
3. Confirm results are identical
4. Measure time savings (~23 min saved)

---

## Phase 2: Stable Trial Period Detection

### Objective
For each neuron, identify a contiguous subset of trials where the neuron's recording is stable, and use only those trials for FR estimation, sparsity calculation, and task modulation analysis.

### Rationale
- A neuron with 10 spikes in trials 1-5, 0 spikes in trials 6-100, and 500 spikes in trials 101-200 has:
  - Current "active period" FR: biased low (first spike to last spike spans trials 1-200)
  - Actual stable FR: much higher if computed only over trials 101-200
- Using the stable period improves:
  - FR threshold accuracy
  - Sparsity estimates
  - Task modulation detection (fewer noisy trials)

### Approach: Sliding Window Presence Analysis

Inspired by the [SpikeInterface presence ratio](https://spikeinterface.readthedocs.io/en/latest/modules/qualitymetrics.html) and [SpikeSift's dynamic segmentation](https://arxiv.org/html/2504.01604v1), we will:

1. **Divide recording into trial-based windows** (e.g., 20-trial sliding windows)
2. **Compute FR in each window** for each neuron
3. **Identify stable plateau** where FR is consistent (low coefficient of variation)
4. **Select the longest stable period** meeting minimum criteria

### Implementation

#### 2.1 Create `detect_stable_trial_period.m`

```matlab
function stability_info = detect_stable_trial_period(spike_times, spike_clusters, ...
    trial_start_times, trial_end_times, cluster_ids, params)
% DETECT_STABLE_TRIAL_PERIOD Identify stable recording period for each neuron.
%
%   This function uses a sliding window approach to detect the longest
%   contiguous period where each neuron's firing rate is stable.
%
%   ALGORITHM:
%   1. For each neuron, compute FR in sliding windows of N trials
%   2. Calculate local coefficient of variation (CV) of FR
%   3. Find longest run of consecutive windows with CV < threshold
%   4. Return trial indices and metrics for the stable period
%
%   INPUTS:
%       spike_times      - [nSpikes × 1] spike times in seconds
%       spike_clusters   - [nSpikes × 1] cluster ID for each spike
%       trial_start_times- [nTrials × 1] start time of each trial
%       trial_end_times  - [nTrials × 1] end time of each trial
%       cluster_ids      - [nClusters × 1] cluster IDs to analyze
%       params           - struct with fields:
%           .window_size     - Number of trials per window (default: 20)
%           .step_size       - Window step in trials (default: 5)
%           .cv_threshold    - Max CV for "stable" (default: 0.5)
%           .min_fr_threshold- Minimum FR to consider present (default: 0.1 sp/s)
%           .min_stable_trials- Minimum trials in stable period (default: 50)
%
%   OUTPUTS:
%       stability_info   - [nClusters × 1] struct array with fields:
%           .cluster_id          - Cluster ID
%           .stable_trial_start  - First trial index of stable period
%           .stable_trial_end    - Last trial index of stable period
%           .n_stable_trials     - Number of trials in stable period
%           .stable_fr           - Mean FR over stable period (sp/s)
%           .stable_sparsity     - Sparsity over stable period
%           .presence_ratio      - Fraction of windows with firing
%           .is_stable           - Boolean: meets stability criteria
%           .stability_score     - Continuous score (0-1, higher = more stable)
%           .window_frs          - [nWindows × 1] FR in each window
%           .window_cvs          - [nWindows × 1] CV in each window
%           .exclusion_reason    - If not stable, why

% Default parameters
if nargin < 6 || isempty(params)
    params = struct();
end
params = set_defaults(params);

nClusters = length(cluster_ids);
nTrials = length(trial_start_times);

% Initialize output
stability_info = struct(...
    'cluster_id', num2cell(cluster_ids), ...
    'stable_trial_start', 1, ...
    'stable_trial_end', nTrials, ...
    'n_stable_trials', nTrials, ...
    'stable_fr', 0, ...
    'stable_sparsity', 0, ...
    'presence_ratio', 0, ...
    'is_stable', false, ...
    'stability_score', 0, ...
    'window_frs', [], ...
    'window_cvs', [], ...
    'exclusion_reason', '');

% Define sliding windows
window_starts = 1:params.step_size:(nTrials - params.window_size + 1);
nWindows = length(window_starts);

for i_cluster = 1:nClusters
    cluster_id = cluster_ids(i_cluster);
    cluster_spike_times = spike_times(spike_clusters == cluster_id);

    if isempty(cluster_spike_times)
        stability_info(i_cluster).exclusion_reason = 'No spikes';
        continue;
    end

    % Compute FR in each sliding window
    window_frs = zeros(nWindows, 1);
    for i_win = 1:nWindows
        win_trials = window_starts(i_win):(window_starts(i_win) + params.window_size - 1);
        win_start = trial_start_times(win_trials(1));
        win_end = trial_end_times(win_trials(end));
        win_duration = win_end - win_start;

        win_spikes = sum(cluster_spike_times >= win_start & cluster_spike_times < win_end);
        window_frs(i_win) = win_spikes / win_duration;
    end

    stability_info(i_cluster).window_frs = window_frs;

    % Calculate presence ratio (fraction of windows with meaningful firing)
    presence_ratio = mean(window_frs >= params.min_fr_threshold);
    stability_info(i_cluster).presence_ratio = presence_ratio;

    % Calculate local CV using rolling window of 3 consecutive windows
    window_cvs = zeros(nWindows - 2, 1);
    for i_win = 1:(nWindows - 2)
        local_frs = window_frs(i_win:(i_win + 2));
        if mean(local_frs) > 0
            window_cvs(i_win) = std(local_frs) / mean(local_frs);
        else
            window_cvs(i_win) = inf;  % No firing = unstable
        end
    end
    stability_info(i_cluster).window_cvs = window_cvs;

    % Find longest stable run (consecutive windows with CV < threshold)
    is_stable_window = window_cvs < params.cv_threshold;
    [longest_run_start, longest_run_length] = find_longest_run(is_stable_window);

    if longest_run_length > 0
        % Convert window indices to trial indices
        stable_win_start = longest_run_start;
        stable_win_end = longest_run_start + longest_run_length - 1;

        stable_trial_start = window_starts(stable_win_start);
        stable_trial_end = min(window_starts(stable_win_end) + params.window_size - 1, nTrials);
        n_stable_trials = stable_trial_end - stable_trial_start + 1;

        % Check minimum trial requirement
        if n_stable_trials >= params.min_stable_trials
            stability_info(i_cluster).stable_trial_start = stable_trial_start;
            stability_info(i_cluster).stable_trial_end = stable_trial_end;
            stability_info(i_cluster).n_stable_trials = n_stable_trials;
            stability_info(i_cluster).is_stable = true;

            % Compute FR and sparsity over stable period only
            stable_start_time = trial_start_times(stable_trial_start);
            stable_end_time = trial_end_times(stable_trial_end);
            stable_spikes = cluster_spike_times(cluster_spike_times >= stable_start_time & ...
                                                 cluster_spike_times < stable_end_time);
            stable_duration = stable_end_time - stable_start_time;

            stability_info(i_cluster).stable_fr = length(stable_spikes) / stable_duration;

            % Sparsity over stable period
            bin_edges = stable_start_time:0.1:stable_end_time;
            if length(bin_edges) > 1
                counts = histcounts(stable_spikes, bin_edges);
                stability_info(i_cluster).stable_sparsity = sum(counts == 0) / length(counts);
            end

            % Stability score (0-1): combines presence ratio and run length
            stability_info(i_cluster).stability_score = ...
                presence_ratio * (n_stable_trials / nTrials);
        else
            stability_info(i_cluster).exclusion_reason = sprintf(...
                'Stable period too short (%d < %d trials)', ...
                n_stable_trials, params.min_stable_trials);
        end
    else
        stability_info(i_cluster).exclusion_reason = 'No stable period found (CV always high)';
    end
end

end

function params = set_defaults(params)
    if ~isfield(params, 'window_size'), params.window_size = 20; end
    if ~isfield(params, 'step_size'), params.step_size = 5; end
    if ~isfield(params, 'cv_threshold'), params.cv_threshold = 0.5; end
    if ~isfield(params, 'min_fr_threshold'), params.min_fr_threshold = 0.1; end
    if ~isfield(params, 'min_stable_trials'), params.min_stable_trials = 50; end
end

function [start_idx, run_length] = find_longest_run(bool_array)
    % Find the longest consecutive run of true values
    if isempty(bool_array) || ~any(bool_array)
        start_idx = 0;
        run_length = 0;
        return;
    end

    % Find run starts and ends
    d = diff([0; bool_array(:); 0]);
    run_starts = find(d == 1);
    run_ends = find(d == -1) - 1;
    run_lengths = run_ends - run_starts + 1;

    [run_length, max_idx] = max(run_lengths);
    start_idx = run_starts(max_idx);
end
```

#### 2.2 Integrate with `apply_neuron_screening.m`

**Add after loading spike data (around line 70):**
```matlab
% Detect stable trial periods for each neuron
fprintf('apply_neuron_screening: Detecting stable trial periods...\n');
stability_params = struct(...
    'window_size', 20, ...
    'step_size', 5, ...
    'cv_threshold', 0.5, ...
    'min_fr_threshold', 0.1, ...
    'min_stable_trials', 50);

% Extract trial timing
trial_start_times = [session_data.trialInfo.trialStartTime]';
trial_end_times = [session_data.trialInfo.trialEndTime]';

stability_info = detect_stable_trial_period(...
    session_data.spikes.times, ...
    session_data.spikes.clusters, ...
    trial_start_times, ...
    trial_end_times, ...
    cluster_ids, ...
    stability_params);

% Store in session_data for caching
session_data.stability_info = stability_info;
```

**Modify FR/sparsity calculation (lines 145-180):**
```matlab
% Use stable period if available, otherwise fall back to active period
if stability_info(i).is_stable
    active_start = trial_start_times(stability_info(i).stable_trial_start);
    active_end = trial_end_times(stability_info(i).stable_trial_end);
    screening_info(i).used_stable_period = true;
    screening_info(i).stable_trial_range = [stability_info(i).stable_trial_start, ...
                                            stability_info(i).stable_trial_end];
else
    % Fall back to first-to-last spike
    active_start = spike_times(1);
    active_end = spike_times(end);
    screening_info(i).used_stable_period = false;
    screening_info(i).stable_trial_range = [];
end
```

#### 2.3 Modify `compute_task_modulation.m`

**Add optional input for stable trial indices:**
```matlab
function modulation_info = compute_task_modulation(session_data, core_data, ...
    condition_defs, stability_info)
% ... existing header ...
%
%   INPUTS:
%       stability_info - (Optional) struct array from detect_stable_trial_period
%                        If provided, only trials in stable period are used

% Handle optional stability_info
use_stable_trials = nargin >= 4 && ~isempty(stability_info);
```

**Modify trial selection (around line 112):**
```matlab
for loc = 1:4
    loc_trials = find(targetLocIdx == loc);

    % Filter to stable trials if available
    if use_stable_trials && stability_info(i_neuron).is_stable
        stable_range = stability_info(i_neuron).stable_trial_start:...
                       stability_info(i_neuron).stable_trial_end;
        loc_trials = intersect(loc_trials, stable_range);
    end

    if length(loc_trials) < 5; continue; end
    % ... rest of existing code ...
end
```

#### 2.4 Cache Stability Info

Add caching similar to `cached_modulation_info`:
```matlab
% Check for cached stability_info
if isfield(session_data, 'stability_info') && ~isempty(session_data.stability_info)
    cached_stability = session_data.stability_info;
    if length(cached_stability) == nNeurons && ...
            isequal([cached_stability.cluster_id], cluster_ids(:)')
        fprintf('apply_neuron_screening: Using cached stability_info.\n');
        stability_info = cached_stability;
    end
end
```

### Testing

#### Unit Tests
1. Create synthetic neuron with known stable period (e.g., silence for trials 1-50, firing for 51-200)
2. Verify `detect_stable_trial_period` correctly identifies trials 51-200
3. Confirm FR computed over stable period matches expected value

#### Integration Tests
1. Run full pipeline on one session
2. Compare screening results with/without stable period detection
3. Manually inspect neurons flagged as "unstable" to validate

#### Validation Metrics
- Count of neurons where stable period differs significantly from active period
- Compare task modulation detection rates with/without stable period filtering
- Visual inspection of stability diagnostic plots

---

## Phase 3: Diagnostic Visualization

### Objective
Create diagnostic plots to visualize neuron stability and validate the detection algorithm.

### Implementation

#### 3.1 Create `plot_neuron_stability.m`

Generate a per-session PDF showing:
1. **Raster plots** for each neuron with stable period highlighted
2. **Sliding window FR** over time with stability threshold
3. **Summary statistics**: presence ratio, stability score, stable vs. total trials

#### 3.2 Add to `run_4factors_analysis.m`

```matlab
% Generate stability diagnostic plots (optional)
if generate_stability_plots
    plot_neuron_stability(session_data, stability_info, output_dir);
end
```

---

## Implementation Order and Timeline

### Phase 1: Friedman Caching (Priority: High, Complexity: Low)
- **Files to modify**: `screen_sc_neurons.m`, `apply_neuron_screening.m`
- **New files**: None
- **Expected time savings**: ~23 min per full screening run
- **Risk**: Low (similar pattern to existing `modulation_info` caching)

### Phase 2: Stable Period Detection (Priority: High, Complexity: Medium)
- **Files to modify**: `apply_neuron_screening.m`, `compute_task_modulation.m`
- **New files**: `detect_stable_trial_period.m`
- **Expected benefit**: More accurate FR/sparsity estimates, better task modulation detection
- **Risk**: Medium (requires validation, parameter tuning)

### Phase 3: Diagnostics (Priority: Medium, Complexity: Low)
- **Files to modify**: `run_4factors_analysis.m`
- **New files**: `plot_neuron_stability.m`
- **Expected benefit**: Easier validation and debugging
- **Risk**: Low

---

## Configuration Parameters

### Recommended Defaults
```matlab
% In define_task_conditions.m, add:
stability_detection = struct(...
    'enabled', true, ...
    'window_size_trials', 20, ...     % ~2-3 min of data
    'step_size_trials', 5, ...        % 75% overlap
    'cv_threshold', 0.5, ...          % Max CV for stable
    'min_fr_threshold', 0.1, ...      % sp/s to count as present
    'min_stable_trials', 50, ...      % Minimum for valid period
    'fallback_to_active', true ...    % Use active period if no stable found
);
```

### Parameter Sensitivity
- **window_size**: Larger = smoother but may miss short stable periods
- **cv_threshold**: Lower = stricter stability requirement
- **min_stable_trials**: Higher = more data required, fewer neurons pass

---

## Success Criteria

1. **Friedman caching**: Second screening run completes in <30% of original time
2. **Stable period detection**:
   - ≥95% of neurons have a detected stable period
   - FR estimates differ by >10% for ≥20% of neurons (indicating meaningful correction)
3. **Task modulation**: Detection rate comparable or improved vs. baseline
4. **No regressions**: Existing screening results are reproducible when stability detection is disabled

---

## References

1. [SpikeInterface Quality Metrics](https://spikeinterface.readthedocs.io/en/latest/modules/qualitymetrics.html) - Presence ratio and other quality metrics
2. [Hill et al. (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3123734/) - Quality metrics for spike sorting
3. [Kilosort4 (Nature Methods, 2024)](https://www.nature.com/articles/s41592-024-02232-7) - Modern drift correction
4. [SpikeSift (2025)](https://arxiv.org/html/2504.01604v1) - Dynamic segmentation for drift resilience
5. [MEDiCINe (eNeuro, 2025)](https://www.eneuro.org/content/12/3/ENEURO.0529-24.2025) - Motion correction for neural recordings
6. [Multi-day Neuron Tracking (eLife, 2024)](https://elifesciences.org/articles/92495) - Tracking across sessions
