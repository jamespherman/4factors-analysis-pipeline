# Neuron Present Trials: Codebase Review

## Executive Summary

This document reviews the codebase to identify where and how **"neuron present trials"** (the subset of trials over which a given neuron/cluster is stably present) should be used. The recent overhaul of neuron screening introduced `detect_stable_trial_period.m`, which identifies stable recording periods for each neuron. However, this vital trial-level selection mask is not yet consistently propagated to downstream analyses, potentially diluting or artificially amplifying metrics like ROC AUC values.

## 1. What Are "Neuron Present Trials"?

### Definition

A **"neuron present trial"** is a trial during which a neuron (Kilosort cluster) is stably present and reliably recorded. This is determined by:

1. **Stable Trial Period Detection** (`detect_stable_trial_period.m`): Uses a sliding window approach to find the longest contiguous period where each neuron's firing rate is stable (CV < threshold across consecutive windows).

2. **NaN-based Trial Exclusion**: Trials where the alignment event time is NaN (missing event) result in NaN firing rate values in `core_data.spikes.<event>.rates`.

### Data Structures

| Field | Location | Description |
|-------|----------|-------------|
| `stability_info(i).stable_trial_start` | `session_data.stability_info` | First trial index of stable period |
| `stability_info(i).stable_trial_end` | `session_data.stability_info` | Last trial index of stable period |
| `stability_info(i).n_stable_trials` | `session_data.stability_info` | Number of trials in stable period |
| `stability_info(i).is_stable` | `session_data.stability_info` | Boolean: neuron has valid stable period |
| `core_data.spikes.<event>.rates` | `session_data.analysis.core_data` | `[n_neurons x n_trials x n_time_bins]` - NaN where event missing |

### Important Note on Trial Indexing

The `stable_trial_start` and `stable_trial_end` indices refer to **ALL trials** in the session (indices into `session_data.eventTimes.trialBegin`), not just the gSac_4factors trials used in `core_data`. A mapping is needed to convert these to the filtered trial indices used in analyses.

## 2. Current Usage of Neuron Present Trials

### 2.1 Screening Functions (Correct Usage)

#### `apply_neuron_screening.m` (lines 176-253)
**Status: CORRECT**

- Uses `stability_info` to compute mean FR and sparsity over stable period
- Stores `used_stable_period` and `stable_trial_range` in `screening_info`
- Falls back to first-to-last spike period if no stable period detected

```matlab
if exist('stability_info', 'var') && stability_info(i).is_stable
    % Use stable period (trial-based timing)
    stable_trial_start = stability_info(i).stable_trial_start;
    stable_trial_end = stability_info(i).stable_trial_end;
    active_start = session_data.eventTimes.trialBegin(stable_trial_start);
    active_end = session_data.eventTimes.trialEnd(stable_trial_end);
```

#### `screen_sc_neurons.m` (lines 210-246)
**Status: PARTIALLY CORRECT**

- Computes epoch firing rates with `valid_trials_mask = ~isnan(epoch_event_times)`
- Removes NaN trials before computing metrics
- **Does NOT use stability_info** - uses all memory-guided saccade trials

#### `compute_task_modulation.m` (lines 133-145)
**Status: PARTIALLY CORRECT**

- Removes NaN trials per neuron/location: `valid_trials = ~any(isnan(neuron_rates), 2)`
- **Does NOT use stability_info**

### 2.2 Analysis Functions (Need Updates)

#### `analyze_window_roc.m` (lines 238-240)
**Status: NEEDS UPDATE**

Current behavior:
- Removes NaN trials after extracting condition data
- Does NOT account for per-neuron stable trial periods

```matlab
cond1_data = neuron_rates(cond1_mask);
cond2_data = neuron_rates(cond2_mask);
cond1_data = cond1_data(~isnan(cond1_data));  % Only NaN removal
cond2_data = cond2_data(~isnan(cond2_data));
```

**Impact**: ROC AUC values may be computed using trials where the neuron wasn't stably recorded, potentially diluting or amplifying true effects.

#### `compute_preferred_locations.m` (lines 101-127)
**Status: NEEDS UPDATE**

Current behavior:
- Uses `'omitnan'` flag when computing mean responses
- Does NOT track per-location trial counts
- Does NOT use stability_info

```matlab
location_responses_visual(:, loc) = mean(visual_window_rates(:, loc_mask), 2, 'omitnan');
```

**Impact**: Preferred location estimates may be biased if computed from trials where the neuron wasn't stably present.

#### `analyze_anova.m` (lines 99-135)
**Status: NEEDS UPDATE**

Current behavior:
- Filters by trial mask (task conditions)
- Filters by valid factor values (non-empty)
- Does NOT filter by per-neuron stable trial periods
- Line 182: Skips neurons with >50% NaN, but this is per-bin not per-neuron-presence

**Impact**: ANOVA results may include data from trials where neurons weren't stably recorded.

#### `analyze_baseline_comparison.m` (lines 110-146)
**Status: NEEDS UPDATE**

Current behavior:
- Selects data based on trial mask (condition)
- No explicit NaN removal before `arrayROC`
- Does NOT use stability_info

```matlab
sub_data = event_data.rates(:, trial_mask, :);
```

**Impact**: Baseline comparison may include trials where neurons weren't stably present.

### 2.3 Export Functions

#### `export_for_rnn.m` (lines 166-191)
**Status: NEEDS UPDATE**

Current behavior:
- Uses all valid gSac_4factors trials for all neurons
- No per-neuron trial filtering based on stable periods

**Impact**: Exported data for RNN training may include trials where specific neurons weren't stably recorded.

## 3. Summary of Files Requiring Updates

| File | Priority | Current Handling | Required Change |
|------|----------|-----------------|-----------------|
| `analyze_window_roc.m` | **HIGH** | NaN removal only | Add per-neuron stable trial filtering |
| `compute_preferred_locations.m` | **HIGH** | `'omitnan'` flag | Add per-neuron stable trial filtering + track counts |
| `compute_task_modulation.m` | **MEDIUM** | NaN removal only | Add per-neuron stable trial filtering |
| `analyze_anova.m` | **MEDIUM** | Trial mask only | Add per-neuron stable trial filtering |
| `analyze_baseline_comparison.m` | **MEDIUM** | No filtering | Add per-neuron stable trial filtering |
| `export_for_rnn.m` | **LOW** | No filtering | Add per-neuron stable trial filtering |
| `screen_sc_neurons.m` | **LOW** | NaN removal | Consider integrating stability_info |

## 4. Data Flow and Trial Index Mapping

### Current Flow
```
Raw Session Data
      |
      v
detect_stable_trial_period() --> stability_info (ALL trial indices)
      |
      v
prepare_core_data() --> core_data.spikes.<event>.rates (gSac_4factors trial indices)
      |
      v
analyze_window_roc() etc. (gSac_4factors trial indices)
```

### Key Challenge: Index Mapping

`stability_info` uses indices into ALL session trials, but `core_data` only includes gSac_4factors trials. A mapping function is needed:

```matlab
% all_trial_indices: indices into session_data (all trials)
% fourfactors_trial_indices: indices used in core_data

% For a given neuron i:
stable_start = stability_info(i).stable_trial_start;  % Index into ALL trials
stable_end = stability_info(i).stable_trial_end;

% Map to core_data indices:
% Trial j in core_data corresponds to fourfactors_trial_indices(j) in session_data
stable_trial_mask = (fourfactors_trial_indices >= stable_start) & ...
                    (fourfactors_trial_indices <= stable_end);
```

## 5. Recommended Implementation Strategy

### 5.1 Create a Utility Function

Create `get_neuron_present_trials.m` that returns a per-neuron trial mask for core_data:

```matlab
function trial_masks = get_neuron_present_trials(stability_info, fourfactors_trial_indices)
% Returns: [n_neurons x n_trials] logical matrix
```

### 5.2 Update Analysis Functions

Each analysis function should:
1. Accept `stability_info` as an optional input
2. Call `get_neuron_present_trials()` to get per-neuron masks
3. Apply the mask before computing metrics
4. Track and report trial counts used

### 5.3 Update Pipeline

`run_4factors_analysis.m` should:
1. Pass `stability_info` to relevant analysis functions
2. Store `fourfactors_trial_indices` in `core_data` for mapping

## 6. Expected Impact

### Window ROC Analysis
- More accurate ROC AUC estimates
- Consistent trial counts across neurons
- Potential changes in significance assessments

### Preferred Location Computation
- More reliable location preference estimates
- Ability to flag neurons with insufficient trials per location

### ANOVA
- More accurate p-values and F-statistics
- Reduced noise from unstable recording periods

### Population Analyses
- Cleaner signal in aggregate metrics
- More interpretable cross-session comparisons

## 7. Validation Strategy

After implementation:
1. Compare old vs new ROC AUC distributions
2. Check for changes in neuron classification rates
3. Verify trial counts are reasonable per neuron
4. Confirm no regressions in existing analyses
