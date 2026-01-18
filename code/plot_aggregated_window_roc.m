%% plot_aggregated_window_roc.m
%
%   Generates publication-quality summary figures for window-based ROC
%   analysis results. Two figures are created:
%
%   Figure 1: Proportion of neurons with significant factor selectivity,
%             organized as a grid (Rows = Factors, Columns = Epochs).
%             Uses bidirectional bars showing cond1 vs cond2 preference.
%
%   Figure 2: Scatter plots comparing factor selectivity across neurons,
%             supporting Hypothesis 2 testing (SNc subregion specialization).
%
% INPUTS:
%   aggregated_data - A struct containing aggregated data for one brain area.
%                     This script requires the `window_roc` field, which
%                     is structured as:
%
%                     .(epoch_name).(factor_name) = [1 x nSessions] struct array
%
%                     Each struct in the array must contain:
%                       - .auc: [nNeurons x 1] vector of ROC AUC values.
%                       - .p: [nNeurons x 1] vector of p-values.
%                       - .ci: [nNeurons x 2] matrix of confidence intervals.
%                       - .n_trials: [nNeurons x 2] matrix of trial counts.
%
%   brain_area_name - A string with the name of the brain area (e.g., 'SC').
%
% Author: Claude Code
% Date: 2026-01-16
%

function plot_aggregated_window_roc(aggregated_data, brain_area_name)

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures', 'window_roc');
addpath(fullfile(project_root, 'code', 'utils'));

if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

%% Input Validation
if ~isfield(aggregated_data, 'window_roc')
    warning('plot_aggregated_window_roc:no_data', ...
        'Window ROC data is missing from the input struct.');
    return;
end

%% Dynamically Discover Grid Layout
epoch_names = fieldnames(aggregated_data.window_roc);
n_epochs = length(epoch_names);

if n_epochs == 0
    warning('plot_aggregated_window_roc:no_epochs', ...
        'No epochs found in window_roc data.');
    return;
end

% Get factor names from the first epoch
factor_names = fieldnames(aggregated_data.window_roc.(epoch_names{1}));
n_factors = length(factor_names);

if n_factors == 0
    warning('plot_aggregated_window_roc:no_factors', ...
        'No factors found in window_roc data.');
    return;
end

%% Define Colors
colors = richColors;
posColor = colors(7,:);   % Blue - for cond1/high preference (AUC > 0.5)
negColor = colors(10,:);  % Orange - for cond2/low preference (AUC < 0.5)
sig_threshold = 0.05;

%% ========================================================================
%  FIGURE 1: Proportion Significant Grid
%  ========================================================================

fig1 = figure('Position', [100, 100, 200 * n_epochs, 200 * n_factors], ...
    'Color', 'w', 'Name', sprintf('Window ROC Summary - %s', brain_area_name));
h_axes = gobjects(n_factors, n_epochs);

% Pre-calculate all proportions for consistent Y-axis scaling
all_props = cell(n_factors, n_epochs);

for i_factor = 1:n_factors
    factor_name = factor_names{i_factor};

    for i_epoch = 1:n_epochs
        epoch_name = epoch_names{i_epoch};

        % Check if this epoch-factor combination exists
        if ~isfield(aggregated_data.window_roc.(epoch_name), factor_name)
            all_props{i_factor, i_epoch} = struct('pos', 0, 'neg', 0, 'n', 0);
            continue;
        end

        session_data = aggregated_data.window_roc.(epoch_name).(factor_name);

        % Concatenate data across sessions
        all_auc = vertcat(session_data.auc);
        all_p = vertcat(session_data.p);

        % Remove NaN values
        valid_mask = ~isnan(all_auc) & ~isnan(all_p);
        all_auc = all_auc(valid_mask);
        all_p = all_p(valid_mask);

        if isempty(all_auc)
            all_props{i_factor, i_epoch} = struct('pos', 0, 'neg', 0, 'n', 0);
            continue;
        end

        % Calculate proportions
        % Positive: significant AND AUC > 0.5 (prefers cond1/high)
        % Negative: significant AND AUC < 0.5 (prefers cond2/low)
        prop_sig_pos = mean(all_p < sig_threshold & all_auc > 0.5);
        prop_sig_neg = mean(all_p < sig_threshold & all_auc < 0.5);

        all_props{i_factor, i_epoch} = struct('pos', prop_sig_pos, ...
            'neg', prop_sig_neg, 'n', length(all_auc));
    end
end

%% Plot Figure 1
for i_factor = 1:n_factors
    factor_name = factor_names{i_factor};

    for i_epoch = 1:n_epochs
        epoch_name = epoch_names{i_epoch};

        h_axes(i_factor, i_epoch) = mySubPlot([n_factors, n_epochs, ...
            (i_factor-1)*n_epochs + i_epoch]);
        hold on;

        props = all_props{i_factor, i_epoch};

        % Bar plot - positive upward, negative downward
        bar(1, props.pos, 0.6, 'FaceColor', posColor, 'EdgeColor', 'none');
        bar(1, -props.neg, 0.6, 'FaceColor', negColor, 'EdgeColor', 'none');

        % Reference line at zero
        line([0.4, 1.6], [0 0], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.5);

        xlim([0.4, 1.6]);
        set(gca, 'XTick', []);
    end
end

%% Figure 1 Formatting
% Standardize Y-axis limits per row
for i_factor = 1:n_factors
    valid_axes_in_row = h_axes(i_factor, isgraphics(h_axes(i_factor, :)));
    if ~isempty(valid_axes_in_row)
        [~, yLims] = outerLims(valid_axes_in_row);
        max_abs_y = max(abs(yLims));
        if max_abs_y == 0; max_abs_y = 0.5; end
        set(valid_axes_in_row, 'YLim', [-max_abs_y, max_abs_y]);
    end

    for i_epoch = 1:n_epochs
        h_ax = h_axes(i_factor, i_epoch);
        if ~isgraphics(h_ax); continue; end

        % Column headers (top row only)
        if i_factor == 1
            title(h_ax, strrep(epoch_names{i_epoch}, '_', ' '), ...
                'FontWeight', 'normal');
        end

        % Row labels (first column only)
        if i_epoch == 1
            ylabel(h_ax, strrep(factor_names{i_factor}, '_', ' '));
        else
            set(h_ax, 'YTickLabel', []);
        end
    end
end

% General formatting
all_valid_axes = h_axes(isgraphics(h_axes));
set(all_valid_axes, 'TickDir', 'Out', 'Color', 'none', ...
    'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

%% Save Figure 1
fig1_filename = fullfile(figures_dir, ...
    sprintf('window_roc_proportion_significant_%s.pdf', brain_area_name));
pdfSave(fig1_filename, fig1.Position(3:4)/72, fig1);

%% ========================================================================
%  FIGURE 2: Scatter Plots for Hypothesis 2 Testing
%  ========================================================================

% Define scatter plot comparisons
% Each row: {x_factor, y_factor, epoch_for_comparison}
scatter_comparisons = {
    'reward', 'salience', 'visual';
    'reward', 'identity', 'visual';
    'probability', 'salience', 'visual';
    'probability', 'identity', 'visual';
    'reward', 'probability', 'visual';    % Within goal-directed
    'salience', 'identity', 'visual'      % Within stimulus-driven
};

n_scatter_plots = size(scatter_comparisons, 1);
n_cols_scatter = 3;
n_rows_scatter = ceil(n_scatter_plots / n_cols_scatter);

fig2 = figure('Position', [100, 100, 300 * n_cols_scatter, 280 * n_rows_scatter], ...
    'Color', 'w', 'Name', sprintf('Window ROC Scatter - %s', brain_area_name));

for i_plot = 1:n_scatter_plots
    x_factor = scatter_comparisons{i_plot, 1};
    y_factor = scatter_comparisons{i_plot, 2};
    epoch_name = scatter_comparisons{i_plot, 3};

    % Check if data exists for this comparison
    if ~isfield(aggregated_data.window_roc, epoch_name)
        continue;
    end
    if ~isfield(aggregated_data.window_roc.(epoch_name), x_factor) || ...
       ~isfield(aggregated_data.window_roc.(epoch_name), y_factor)
        continue;
    end

    ax = mySubPlot([n_rows_scatter, n_cols_scatter, i_plot]);
    hold on;

    x_data_all = aggregated_data.window_roc.(epoch_name).(x_factor);
    y_data_all = aggregated_data.window_roc.(epoch_name).(y_factor);

    % Concatenate across sessions
    x_auc = vertcat(x_data_all.auc);
    y_auc = vertcat(y_data_all.auc);

    % Remove NaN values (must be valid in both)
    valid_mask = ~isnan(x_auc) & ~isnan(y_auc);
    x_auc = x_auc(valid_mask);
    y_auc = y_auc(valid_mask);

    if isempty(x_auc)
        title(ax, sprintf('%s vs %s (No Data)', x_factor, y_factor));
        xlim([0 1]);
        ylim([0 1]);
        continue;
    end

    % Plot scatter
    scatter(x_auc, y_auc, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.4);

    % Reference lines at 0.5 (chance level)
    line([0.5 0.5], [0 1], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
    line([0 1], [0.5 0.5], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');

    % Unity line
    line([0 1], [0 1], 'Color', [0.85 0.85 0.85], 'LineStyle', '-');

    xlim([0 1]);
    ylim([0 1]);
    axis square;

    xlabel(sprintf('%s AUC', strrep(x_factor, '_', ' ')));
    ylabel(sprintf('%s AUC', strrep(y_factor, '_', ' ')));

    % Add neuron count to title
    title(sprintf('n = %d', length(x_auc)), 'FontWeight', 'normal');
end

% General formatting
set(findall(fig2, 'Type', 'Axes'), 'TickDir', 'Out', 'LineWidth', 1, ...
    'XColor', 'k', 'YColor', 'k');

%% Save Figure 2
fig2_filename = fullfile(figures_dir, ...
    sprintf('window_roc_scatter_%s.pdf', brain_area_name));
pdfSave(fig2_filename, fig2.Position(3:4)/72, fig2);

end
