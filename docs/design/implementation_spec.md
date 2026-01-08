# Implementation Specification Document: 4factors Analysis Pipeline

This document contains technical implementation details for the 4factors analysis pipeline, including data structures, function signatures, and code patterns. It complements the Analysis Rationale Document, which focuses on conceptual design.

---

## Table of Contents

1. [Data Structures](#1-data-structures)
2. [Function Specifications](#2-function-specifications)
3. [Visualization Implementation](#3-visualization-implementation)
4. [File Naming Conventions](#4-file-naming-conventions)

---

## 1. Data Structures

### 1.1 Neuron Metrics Table

Stored in `aggregated_analysis_data.mat`, this table contains per-neuron summary statistics for population-level visualizations.

```matlab
neuron_metrics_table = table(...
    session_id, ...          % e.g., 'Feynman_08_15_2025_SC'
    cluster_id, ...          % e.g., 42
    brain_area, ...          % 'SC' or 'SNc'
    snc_subregion, ...       % 'rvmSNc', 'cdlSNc', or '' for SC
    neuron_type, ...         % 'visually_excited', 'visually_suppressed', etc.
    baseline_fr, ...         % Baseline firing rate (Hz)
    waveform_duration, ...   % Peak-trough duration (ms)
    % Visual epoch ROC values
    visual_reward_roc, ...
    visual_salience_roc, ...
    visual_probability_roc, ...
    visual_identity_roc, ...
    % Saccade epoch ROC values
    saccade_reward_roc, ...
    saccade_salience_roc, ...
    saccade_probability_roc, ...
    saccade_identity_roc, ...
    % Delay epoch ROC values
    delay_reward_roc, ...
    delay_salience_roc, ...
    delay_probability_roc, ...
    delay_identity_roc, ...
    % Post-reward epoch ROC values
    reward_reward_roc, ...
    reward_salience_roc, ...
    reward_probability_roc, ...
    reward_identity_roc ...
);
```

### 1.2 Window-Based ROC Output Structure

```matlab
window_roc.visual.reward.auc      % [n_neurons × 1]
window_roc.visual.reward.p        % [n_neurons × 1]
window_roc.visual.reward.ci       % [n_neurons × 2]
window_roc.visual.salience.auc
window_roc.visual.salience.p
window_roc.visual.salience.ci
% ... etc for all epoch × factor combinations
```

---

## 2. Function Specifications

### 2.1 Window-Based ROC Analysis

```matlab
function window_roc = analyze_window_roc(session_data, conditions, ...
    core_data, config)
% ANALYZE_WINDOW_ROC Compute ROC AUC in fixed time windows.
%
% INPUTS:
%   session_data - Standard session data structure
%   conditions   - Output of define_task_conditions
%   core_data    - Binned spike data from prepare_core_data
%   config       - Configuration struct with fields:
%       .windows: struct array defining epochs
%           .name: 'visual', 'delay', 'saccade', 'reward'
%           .event: 'targetOn', 'fixOff', 'saccadeOnset', 'reward'
%           .window_sc: [start_ms, end_ms] for SC neurons
%           .window_snc: [start_ms, end_ms] for SNc neurons
%       .comparisons: struct array defining factor comparisons
%           .name: 'reward', 'salience', 'probability', 'identity'
%           .cond1: condition mask name for high/condition 1
%           .cond2: condition mask name for low/condition 2
%           .trial_mask: additional trial filter (e.g., 'is_contralateral_target')
%       .n_permutations: number of permutations for significance (default: 1000)
%
% OUTPUT:
%   window_roc: struct with fields for each epoch.factor combination
%       .auc: [n_neurons × 1] ROC AUC values
%       .p: [n_neurons × 1] p-values from permutation test
%       .ci: [n_neurons × 2] 95% confidence intervals
```

### 2.2 Neuron Classification

```matlab
function [neuron_type, classification_stats] = classify_sc_neurons(...
    epoch_frs, trial_masks, p_friedman)
% CLASSIFY_SC_NEURONS Categorize SC neurons by response profile.
%
% INPUTS:
%   epoch_frs: [n_neurons × 4 × n_trials] firing rates
%              Epochs: baseline, visual, delay, saccade
%   trial_masks: struct with fields
%       .is_contralateral: logical mask for contralateral target trials
%   p_friedman: [n_neurons × 1] p-values from Friedman test
%
% OUTPUT:
%   neuron_type: categorical array with values:
%       'visually_excited' - Significant increase in visual epoch (contra)
%       'visually_suppressed' - Significant decrease in visual epoch (contra)
%       'delay_modulated' - Significant change in delay epoch
%       'saccade_modulated' - Significant change in saccade epoch
%       'unmodulated' - No significant modulation
%   classification_stats: struct with detailed statistics per neuron
%
% NOTE: Classification uses Bonferroni-corrected Wilcoxon tests (alpha = 0.05/3)
%       Suppression criterion: significant DECREASE for CONTRALATERAL targets
```

### 2.3 SNc Subregion Parsing

```matlab
function subregion = parse_snc_subregion(grid_hole_str, brain_area)
% PARSE_SNC_SUBREGION Determine SNc subregion from grid coordinates.
%
% INPUT:
%   grid_hole_str: String in format '(X, Y)' from manifest
%   brain_area: String indicating brain area ('SC' or 'SNc')
%
% OUTPUT:
%   subregion: 'rvmSNc', 'cdlSNc', or '' (empty for SC)
%
% RULE (only applied when brain_area == 'SNc'):
%   Y ∈ {3, 4} → 'rvmSNc'
%   Y ∈ {1, 2} → 'cdlSNc'

if ~strcmp(brain_area, 'SNc')
    subregion = '';
    return;
end

coords = sscanf(grid_hole_str, '(%f, %f)');
if length(coords) >= 2
    grid_y = coords(2);
    if grid_y >= 3
        subregion = 'rvmSNc';
    else
        subregion = 'cdlSNc';
    end
else
    subregion = 'unknown';
    warning('Could not parse grid_hole: %s', grid_hole_str);
end
```

---

## 3. Visualization Implementation

### 3.1 Scatter Plot with Hover Tooltips

```matlab
function fig = create_roc_scatter(neuron_table, x_field, y_field, ...
    color_field, title_str)
% CREATE_ROC_SCATTER Create scatter plot with hover tooltips for neuron ID.
%
% INPUTS:
%   neuron_table: table with neuron metrics
%   x_field: string, column name for x-axis values
%   y_field: string, column name for y-axis values
%   color_field: string, column name for color grouping
%   title_str: string, figure title
%
% OUTPUT:
%   fig: figure handle

fig = figure('Color', 'w');

% Extract data
x_data = neuron_table.(x_field);
y_data = neuron_table.(y_field);
groups = neuron_table.(color_field);

% Create scatter
gscatter(x_data, y_data, groups);

% Add reference lines
hold on;
plot([0.5 0.5], ylim, 'k:', 'LineWidth', 0.5);
plot(xlim, [0.5 0.5], 'k:', 'LineWidth', 0.5);

% Configure data cursor
dcm = datacursormode(fig);
dcm.Enable = 'on';
dcm.UpdateFcn = @(~, event) tooltip_callback(event, neuron_table, x_field, y_field);

% Store table in figure
fig.UserData.neuron_table = neuron_table;
fig.UserData.x_field = x_field;
fig.UserData.y_field = y_field;

title(title_str, 'Interpreter', 'none');
xlabel(strrep(x_field, '_', ' '));
ylabel(strrep(y_field, '_', ' '));

end

function txt = tooltip_callback(event, neuron_table, x_field, y_field)
    x = event.Position(1);
    y = event.Position(2);
    
    % Find matching neuron
    tol = 0.005;
    x_vals = neuron_table.(x_field);
    y_vals = neuron_table.(y_field);
    match_idx = find(abs(x_vals - x) < tol & abs(y_vals - y) < tol, 1);
    
    if ~isempty(match_idx)
        txt = sprintf(['Session: %s\n' ...
                       'Cluster: %d\n' ...
                       'Type: %s\n' ...
                       'Subregion: %s\n' ...
                       '%s: %.3f\n' ...
                       '%s: %.3f'], ...
            neuron_table.session_id{match_idx}, ...
            neuron_table.cluster_id(match_idx), ...
            neuron_table.neuron_type{match_idx}, ...
            neuron_table.snc_subregion{match_idx}, ...
            strrep(x_field, '_', ' '), x, ...
            strrep(y_field, '_', ' '), y);
    else
        txt = sprintf('X: %.3f\nY: %.3f', x, y);
    end
end
```

### 3.2 CSV Backup for Scatter Plots

When generating any scatter plot, also export a CSV lookup file:

```matlab
% After creating scatter plot
csv_filename = fullfile(figures_dir, 'scatter_data', ...
    sprintf('%s_neuron_lookup.csv', plot_name));
writetable(neuron_table(:, {'session_id', 'cluster_id', 'brain_area', ...
    'snc_subregion', 'neuron_type', x_field, y_field}), csv_filename);
```

### 3.3 Linking Simultaneous Sessions

```matlab
function paired_sessions = find_simultaneous_sessions(manifest)
% FIND_SIMULTANEOUS_SESSIONS Identify SC-SNc pairs recorded simultaneously.
%
% INPUT:
%   manifest: table from session_manifest.csv
%
% OUTPUT:
%   paired_sessions: struct array with fields:
%       .group_id: session_group_id
%       .sc_session: unique_id for SC session
%       .snc_session: unique_id for SNc session
%       .snc_subregion: 'rvmSNc' or 'cdlSNc'

unique_groups = unique(manifest.session_group_id);
paired_sessions = struct('group_id', {}, 'sc_session', {}, ...
    'snc_session', {}, 'snc_subregion', {});

for i = 1:length(unique_groups)
    group_id = unique_groups{i};
    group_rows = strcmp(manifest.session_group_id, group_id);
    
    if sum(group_rows) == 2  % Paired SC + SNc
        sessions = manifest(group_rows, :);
        sc_idx = strcmp(sessions.brain_area, 'SC');
        snc_idx = strcmp(sessions.brain_area, 'SNc');
        
        if any(sc_idx) && any(snc_idx)
            entry.group_id = group_id;
            entry.sc_session = sessions.unique_id{sc_idx};
            entry.snc_session = sessions.unique_id{snc_idx};
            entry.snc_subregion = parse_snc_subregion(...
                sessions.grid_hole{snc_idx}, 'SNc');
            paired_sessions(end+1) = entry;
        end
    end
end
```

---

## 4. File Naming Conventions

### 4.1 Output Files

| File Type | Location | Naming Pattern |
|-----------|----------|----------------|
| Session analysis | `output/` | `{unique_id}_session_data.mat` |
| Aggregated data | `data/processed/` | `aggregated_analysis_data.mat` |
| Neuron metrics | `data/processed/` | `aggregated_neuron_metrics.mat` |
| Neuron PDFs | `figures/neuron_summaries/` | `{unique_id}_{cluster_id:03d}_summary.pdf` |
| Population figures | `figures/` | `aggregated_{analysis_type}_{brain_area}.pdf` |
| Scatter data | `figures/scatter_data/` | `{plot_name}_neuron_lookup.csv` |

### 4.2 Figure Subdirectories

```
figures/
├── neuron_summaries/     # Per-neuron diagnostic PDFs
├── decoding/             # Population decoding results
├── scatter_data/         # CSV backups for scatter plots
└── screening/            # Neuron screening diagnostic plots
```

---

## Document History

| Date | Author | Description |
|------|--------|-------------|
| 2025-01-07 | Claude (with James Herman) | Initial creation, split from Analysis Rationale Document |
