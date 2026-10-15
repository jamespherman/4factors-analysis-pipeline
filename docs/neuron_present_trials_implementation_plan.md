# Implementation Plan: Neuron Present Trials

## Overview

This document outlines the implementation plan for propagating "neuron present trials" (stable trial periods) throughout the analysis pipeline. The goal is to ensure that all per-neuron analyses only use trials where that neuron was stably recorded.

## Phase 1: Infrastructure

### Task 1.1: Create Utility Function `get_neuron_present_trials.m`

**Location**: `code/utils/get_neuron_present_trials.m`

**Purpose**: Returns a per-neuron trial mask indicating which core_data trials fall within each neuron's stable recording period.

**Interface**:
```matlab
function neuron_trial_masks = get_neuron_present_trials(stability_info, ...
    fourfactors_trial_indices, n_neurons)
%
% INPUTS:
%   stability_info           - Struct array from detect_stable_trial_period
%   fourfactors_trial_indices - Trial indices used in core_data (into session_data)
%   n_neurons                - Number of neurons
%
% OUTPUT:
%   neuron_trial_masks       - [n_neurons x n_trials] logical matrix
%                              True where trial is within neuron's stable period
```

**Logic**:
1. For each neuron, get `stable_trial_start` and `stable_trial_end`
2. Map these to the filtered gSac_4factors trial indices
3. Return logical mask for each neuron

### Task 1.2: Store `fourfactors_trial_indices` in `core_data`

**File**: `code/prepare_core_data.m`

**Change**: Add the trial indices to the output struct so downstream functions can perform the mapping.

```matlab
core_data.trial_indices = fourfactors_trial_indices;
```

## Phase 2: High-Priority Analysis Updates

### Task 2.1: Update `analyze_window_roc.m`

**Purpose**: Use neuron present trials when computing ROC AUC values.

**Changes**:
1. Accept optional `stability_info` parameter
2. Compute per-neuron trial masks
3. Apply mask in the per-neuron loop (lines 197-265)
4. Track number of trials used per neuron

**Modified Loop (pseudocode)**:
```matlab
for i_neuron = 1:n_neurons
    % Get stable trial mask for this neuron
    if use_stable_trials
        stable_mask = neuron_trial_masks(i_neuron, :)';
    else
        stable_mask = true(n_trials, 1);
    end

    % Apply stable mask to condition masks
    cond1_mask_stable = cond1_mask & stable_mask;
    cond2_mask_stable = cond2_mask & stable_mask;

    % Extract data using stable masks
    cond1_data = neuron_rates(cond1_mask_stable);
    cond2_data = neuron_rates(cond2_mask_stable);

    % Remove NaN (for missing event times)
    cond1_data = cond1_data(~isnan(cond1_data));
    cond2_data = cond2_data(~isnan(cond2_data));

    % ... rest of ROC computation
end
```

### Task 2.2: Update `compute_preferred_locations.m`

**Purpose**: Use neuron present trials when computing preferred locations.

**Changes**:
1. Accept optional `stability_info` and `trial_indices` parameters
2. For each neuron, only use trials within stable period
3. Track trial counts per location per neuron
4. Add minimum trial threshold per location

**Key modification** (per-neuron per-location):
```matlab
for i_neuron = 1:n_neurons
    for loc = 1:4
        loc_mask = (targetLocIdx == loc);
        if use_stable_trials
            loc_mask = loc_mask & neuron_trial_masks(i_neuron, :)';
        end
        % Also remove NaN trials
        neuron_rates = visual_window_rates(i_neuron, :);
        valid_mask = loc_mask & ~isnan(neuron_rates');

        n_valid = sum(valid_mask);
        if n_valid >= min_trials_per_loc
            location_responses_visual(i_neuron, loc) = mean(neuron_rates(valid_mask));
            trial_counts(i_neuron, loc) = n_valid;
        end
    end
end
```

**Output additions**:
- `preferred_locations.trial_counts` - [n_neurons x 4] trial counts per location

## Phase 3: Medium-Priority Analysis Updates

### Task 3.1: Update `compute_task_modulation.m`

**Purpose**: Use neuron present trials when computing task modulation.

**Changes**:
1. Accept optional `stability_info` and `trial_indices` parameters
2. Apply per-neuron stable trial mask before computing modulation
3. Update trial count tracking

### Task 3.2: Update `analyze_anova.m`

**Purpose**: Use neuron present trials in ANOVA analysis.

**Changes**:
1. Accept optional `stability_info` parameter
2. For each neuron, filter trials to stable period before ANOVA
3. Skip neurons with insufficient stable trials

**Key modification**:
```matlab
for i_neuron = 1:n_neurons
    for i_bin = 1:n_bins
        y = squeeze(fr_filt(i_neuron, :, i_bin))';

        % Apply neuron-specific stable trial mask
        if use_stable_trials
            neuron_stable = neuron_trial_masks(i_neuron, valid_trials_mask);
            y_stable = y(neuron_stable);
            group_vars_stable = cellfun(@(g) g(neuron_stable), group_vars, 'UniformOutput', false);
        else
            y_stable = y;
            group_vars_stable = group_vars;
        end

        % ... ANOVA on y_stable with group_vars_stable
    end
end
```

### Task 3.3: Update `analyze_baseline_comparison.m`

**Purpose**: Use neuron present trials in baseline comparisons.

**Changes**:
1. Accept optional `stability_info` parameter
2. Apply per-neuron stable trial mask before arrayROC
3. Handle cases where insufficient trials remain

## Phase 4: Pipeline Integration

### Task 4.1: Update `run_4factors_analysis.m`

**Changes**:
1. After screening, store `fourfactors_trial_indices` in `core_data` (via prepare_core_data)
2. Pass `stability_info` to:
   - `compute_preferred_locations()` (line 637)
   - `analyze_window_roc()` (line 654)
   - `analyze_baseline_comparison()` (line 545)
   - `analyze_anova()` (line 702)
3. Force re-run of analyses when screening changes

### Task 4.2: Update Analysis Re-run Logic

When `stability_info` changes or screening is re-run, analyses that depend on neuron present trials should be invalidated and re-run. Add a version/hash mechanism or simply set `force_rerun.analyses = true` after screening.

## Phase 5: Low-Priority Updates (Optional)

### Task 5.1: Update `export_for_rnn.m`

**Purpose**: Include neuron present trial information in RNN exports.

**Options**:
1. Export only trials within stable period per neuron (results in per-neuron trial subsets)
2. Export a trial mask per neuron
3. Mark trials outside stable period with a flag

**Recommended approach**: Option 2 - export `neuron_stable_trial_mask` [n_neurons x n_trials] logical matrix.

### Task 5.2: Update `screen_sc_neurons.m`

**Purpose**: Optionally use stability_info for epoch FR calculations.

**Note**: Lower priority since this function computes its own per-location trial filtering. May not significantly impact results.

## Implementation Order

1. **Task 1.1**: Create `get_neuron_present_trials.m` utility
2. **Task 1.2**: Store trial indices in `core_data`
3. **Task 2.1**: Update `analyze_window_roc.m`
4. **Task 2.2**: Update `compute_preferred_locations.m`
5. **Task 4.1**: Update `run_4factors_analysis.m` to pass stability_info
6. **Task 3.1**: Update `compute_task_modulation.m`
7. **Task 3.2**: Update `analyze_anova.m`
8. **Task 3.3**: Update `analyze_baseline_comparison.m`
9. **Task 5.1**: Update `export_for_rnn.m` (optional)

## Testing Strategy

### Unit Tests
1. `get_neuron_present_trials.m`: Test mapping logic with known indices
2. Each updated function: Verify trial counts match expectations

### Integration Tests
1. Run full pipeline on one SC and one SNc session
2. Compare old vs new results:
   - ROC AUC distributions
   - Preferred location assignments
   - Number of significant neurons

### Validation Metrics
1. Per-neuron trial counts should be <= stable trial counts
2. No negative trial counts
3. Sum of per-condition trials should match total trials within stable period

## Rollback Plan

If issues are discovered:
1. All changes are backward compatible (new parameters are optional)
2. Set `use_stable_trials = false` to revert to old behavior
3. Re-run pipeline with `force_rerun.analyses = true`

## Expected Timeline

- Phase 1 (Infrastructure): Core utility functions
- Phase 2 (High Priority): Window ROC and preferred locations
- Phase 3 (Medium Priority): ANOVA and baseline comparison
- Phase 4 (Integration): Pipeline updates
- Phase 5 (Optional): Export and SC screening
