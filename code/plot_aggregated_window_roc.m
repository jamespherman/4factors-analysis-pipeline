%% plot_aggregated_window_roc.m
%
%   Generates publication-quality summary figures for window-based ROC
%   analysis results. Two figures are created:
%
%   Figure 1: Proportion of neurons with significant factor selectivity,
%             organized as a grid (Rows = Factors, Columns = Epochs).
%             Uses bidirectional bars showing cond1 vs cond2 preference.
%             For SC data, shows separate bars for 'classic' and 'interneuron'.
%
%   Figure 2: Scatter plots comparing factor selectivity across neurons,
%             supporting Hypothesis 2 testing (SNc subregion specialization).
%             Points are color-coded by significance category (neither, x-only,
%             y-only, or both factors significant based on 95% CI).
%             For SC data, uses circles for 'classic' and stars for 'interneuron'.
%             Each panel has its own legend with counts specific to that comparison.
%
%   Note: Neuron counts may vary slightly between panels because some neurons
%   have NaN values for certain factors (due to insufficient trial counts in
%   analyze_window_roc.m). Each panel shows only neurons with valid data for
%   BOTH factors being compared.
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
%                     For SC data, also requires `neuron_metrics_table` with:
%                       - neuron_type: 'classic', 'interneuron', or 'excluded'
%                       - is_selected: logical (task-modulated neurons)
%                       - auc_{epoch}_{factor}: ROC AUC values
%                       - ci_lower_{epoch}_{factor}, ci_upper_{epoch}_{factor}: CI bounds
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

% Significance category colors (following example image style)
color_neither = [0.6 0.6 0.6];      % Gray for neither significant
color_x_only = colors(7,:);          % Blue for x-factor only
color_y_only = colors(10,:);         % Orange/red for y-factor only
color_both = colors(4,:);            % Magenta/purple for both

sig_threshold = 0.05;

%% Determine if we should differentiate SC neuron types
use_sc_neuron_types = strcmpi(brain_area_name, 'SC') && ...
    isfield(aggregated_data, 'neuron_metrics_table') && ...
    ~isempty(aggregated_data.neuron_metrics_table) && ...
    ismember('neuron_type', aggregated_data.neuron_metrics_table.Properties.VariableNames);

% Pre-extract SC neuron type data if available
if use_sc_neuron_types
    nmt = aggregated_data.neuron_metrics_table;
    % Filter for selected neurons (classic + interneuron, not excluded)
    selected_mask = nmt.is_selected;
    nmt_selected = nmt(selected_mask, :);

    % Identify classic vs interneuron
    is_classic = strcmp(nmt_selected.neuron_type, 'classic');
    is_interneuron = strcmp(nmt_selected.neuron_type, 'interneuron');

    fprintf('SC neuron types (total selected): %d classic, %d interneuron\n', ...
        sum(is_classic), sum(is_interneuron));
end

%% ========================================================================
%  FIGURE 1: Proportion Significant Grid
%  ========================================================================

if use_sc_neuron_types
    % For SC: wider figure to accommodate side-by-side bars
    fig1 = figure('Position', [100, 100, 280 * n_epochs, 200 * n_factors], ...
        'Color', 'w', 'Name', sprintf('Window ROC Summary - %s', brain_area_name));
else
    fig1 = figure('Position', [100, 100, 200 * n_epochs, 200 * n_factors], ...
        'Color', 'w', 'Name', sprintf('Window ROC Summary - %s', brain_area_name));
end
h_axes = gobjects(n_factors, n_epochs);

% Pre-calculate all proportions for consistent Y-axis scaling
all_props = cell(n_factors, n_epochs);

for i_factor = 1:n_factors
    factor_name = factor_names{i_factor};

    for i_epoch = 1:n_epochs
        epoch_name = epoch_names{i_epoch};

        if use_sc_neuron_types
            % Calculate proportions separately for classic and interneuron
            auc_col = sprintf('auc_%s_%s', epoch_name, factor_name);
            p_col = sprintf('p_%s_%s', epoch_name, factor_name);

            if ~ismember(auc_col, nmt_selected.Properties.VariableNames) || ...
               ~ismember(p_col, nmt_selected.Properties.VariableNames)
                all_props{i_factor, i_epoch} = struct(...
                    'classic_pos', 0, 'classic_neg', 0, 'n_classic', 0, ...
                    'intern_pos', 0, 'intern_neg', 0, 'n_intern', 0);
                continue;
            end

            auc_vals = nmt_selected.(auc_col);
            p_vals = nmt_selected.(p_col);

            % Classic neurons
            classic_mask = is_classic & ~isnan(auc_vals) & ~isnan(p_vals);
            classic_auc = auc_vals(classic_mask);
            classic_p = p_vals(classic_mask);

            if ~isempty(classic_auc)
                classic_pos = mean(classic_p < sig_threshold & classic_auc > 0.5);
                classic_neg = mean(classic_p < sig_threshold & classic_auc < 0.5);
            else
                classic_pos = 0; classic_neg = 0;
            end

            % Interneurons
            intern_mask = is_interneuron & ~isnan(auc_vals) & ~isnan(p_vals);
            intern_auc = auc_vals(intern_mask);
            intern_p = p_vals(intern_mask);

            if ~isempty(intern_auc)
                intern_pos = mean(intern_p < sig_threshold & intern_auc > 0.5);
                intern_neg = mean(intern_p < sig_threshold & intern_auc < 0.5);
            else
                intern_pos = 0; intern_neg = 0;
            end

            all_props{i_factor, i_epoch} = struct(...
                'classic_pos', classic_pos, 'classic_neg', classic_neg, ...
                'n_classic', sum(classic_mask), ...
                'intern_pos', intern_pos, 'intern_neg', intern_neg, ...
                'n_intern', sum(intern_mask));
        else
            % Original behavior for non-SC data
            if ~isfield(aggregated_data.window_roc.(epoch_name), factor_name)
                all_props{i_factor, i_epoch} = struct('pos', 0, 'neg', 0, 'n', 0);
                continue;
            end

            session_data = aggregated_data.window_roc.(epoch_name).(factor_name);
            all_auc = vertcat(session_data.auc);
            all_p = vertcat(session_data.p);

            valid_mask = ~isnan(all_auc) & ~isnan(all_p);
            all_auc = all_auc(valid_mask);
            all_p = all_p(valid_mask);

            if isempty(all_auc)
                all_props{i_factor, i_epoch} = struct('pos', 0, 'neg', 0, 'n', 0);
                continue;
            end

            prop_sig_pos = mean(all_p < sig_threshold & all_auc > 0.5);
            prop_sig_neg = mean(all_p < sig_threshold & all_auc < 0.5);

            all_props{i_factor, i_epoch} = struct('pos', prop_sig_pos, ...
                'neg', prop_sig_neg, 'n', length(all_auc));
        end
    end
end

%% Plot Figure 1
% Define bar colors for classic vs interneuron
classic_pos_color = posColor;
classic_neg_color = negColor;
intern_pos_color = posColor * 0.6;   % Darker blue for interneuron
intern_neg_color = negColor * 0.6;   % Darker orange for interneuron

for i_factor = 1:n_factors
    factor_name = factor_names{i_factor};

    for i_epoch = 1:n_epochs
        epoch_name = epoch_names{i_epoch};

        h_axes(i_factor, i_epoch) = mySubPlot([n_factors, n_epochs, ...
            (i_factor-1)*n_epochs + i_epoch]);
        hold on;

        props = all_props{i_factor, i_epoch};

        if use_sc_neuron_types
            % Side-by-side bars for classic (left) and interneuron (right)
            bar_width = 0.35;

            % Classic bars (position 0.8)
            bar(0.8, props.classic_pos, bar_width, 'FaceColor', classic_pos_color, 'EdgeColor', 'none');
            bar(0.8, -props.classic_neg, bar_width, 'FaceColor', classic_neg_color, 'EdgeColor', 'none');

            % Interneuron bars (position 1.2)
            bar(1.2, props.intern_pos, bar_width, 'FaceColor', intern_pos_color, 'EdgeColor', 'none');
            bar(1.2, -props.intern_neg, bar_width, 'FaceColor', intern_neg_color, 'EdgeColor', 'none');

            % Reference line at zero
            line([0.4, 1.6], [0 0], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.5);

            xlim([0.4, 1.6]);
            set(gca, 'XTick', [0.8, 1.2], 'XTickLabel', {'C', 'I'}, 'FontSize', 7);
        else
            % Original single bar plot
            bar(1, props.pos, 0.6, 'FaceColor', posColor, 'EdgeColor', 'none');
            bar(1, -props.neg, 0.6, 'FaceColor', negColor, 'EdgeColor', 'none');

            line([0.4, 1.6], [0 0], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 0.5);

            xlim([0.4, 1.6]);
            set(gca, 'XTick', []);
        end
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

% Add legend for SC neuron types
if use_sc_neuron_types
    % Create legend in first subplot
    axes(h_axes(1, 1));
    hold on;
    h_leg1 = bar(NaN, NaN, 'FaceColor', classic_pos_color, 'EdgeColor', 'none');
    h_leg2 = bar(NaN, NaN, 'FaceColor', intern_pos_color, 'EdgeColor', 'none');
    legend([h_leg1, h_leg2], {'Classic', 'Interneuron'}, ...
        'Location', 'best', 'FontSize', 7, 'Box', 'off');
end

%% Save Figure 1
try
    theme('light');
catch
end

fig1_filename = fullfile(figures_dir, ...
    sprintf('window_roc_proportion_significant_%s.pdf', brain_area_name));
pdfSave(fig1_filename, fig1.Position(3:4)/72, fig1);

%% ========================================================================
%  FIGURE 2: Scatter Plots for Hypothesis 2 Testing
%  ========================================================================

% Define scatter plot comparisons
scatter_comparisons = {
    'reward', 'salience', 'visual';
    'reward', 'identity', 'visual';
    'probability', 'salience', 'visual';
    'probability', 'identity', 'visual';
    'reward', 'probability', 'visual';
    'salience', 'identity', 'visual'
};

n_scatter_plots = size(scatter_comparisons, 1);
n_cols_scatter = 3;
n_rows_scatter = ceil(n_scatter_plots / n_cols_scatter);

fig2 = figure('Position', [100, 100, 350 * n_cols_scatter, 320 * n_rows_scatter], ...
    'Color', 'w', 'Name', sprintf('Window ROC Scatter - %s', brain_area_name));

% Marker settings
marker_size = 40;
legend_marker_size = 25;
axes_bg_color = [1 1 1];  % White background - used for MarkerEdgeColor

for i_plot = 1:n_scatter_plots
    x_factor = scatter_comparisons{i_plot, 1};
    y_factor = scatter_comparisons{i_plot, 2};
    epoch_name = scatter_comparisons{i_plot, 3};

    ax = mySubPlot([n_rows_scatter, n_cols_scatter, i_plot]);
    hold on;

    if use_sc_neuron_types
        % Use neuron_metrics_table for SC data with CI-based significance
        x_auc_col = sprintf('auc_%s_%s', epoch_name, x_factor);
        y_auc_col = sprintf('auc_%s_%s', epoch_name, y_factor);
        x_ci_lower_col = sprintf('ci_lower_%s_%s', epoch_name, x_factor);
        x_ci_upper_col = sprintf('ci_upper_%s_%s', epoch_name, x_factor);
        y_ci_lower_col = sprintf('ci_lower_%s_%s', epoch_name, y_factor);
        y_ci_upper_col = sprintf('ci_upper_%s_%s', epoch_name, y_factor);

        % Check if columns exist
        required_cols = {x_auc_col, y_auc_col, x_ci_lower_col, x_ci_upper_col, ...
                        y_ci_lower_col, y_ci_upper_col};
        cols_exist = all(ismember(required_cols, nmt_selected.Properties.VariableNames));

        if ~cols_exist
            title(ax, sprintf('%s vs %s (No CI Data)', x_factor, y_factor));
            xlim([0 1]); ylim([0 1]);
            continue;
        end

        % Extract data
        x_auc = nmt_selected.(x_auc_col);
        y_auc = nmt_selected.(y_auc_col);
        x_ci_lower = nmt_selected.(x_ci_lower_col);
        x_ci_upper = nmt_selected.(x_ci_upper_col);
        y_ci_lower = nmt_selected.(y_ci_lower_col);
        y_ci_upper = nmt_selected.(y_ci_upper_col);

        % Valid mask (non-NaN in all required fields)
        valid_mask = ~isnan(x_auc) & ~isnan(y_auc) & ...
                     ~isnan(x_ci_lower) & ~isnan(x_ci_upper) & ...
                     ~isnan(y_ci_lower) & ~isnan(y_ci_upper);

        % Determine significance based on CI not containing 0.5
        x_sig = (x_ci_lower > 0.5) | (x_ci_upper < 0.5);
        y_sig = (y_ci_lower > 0.5) | (y_ci_upper < 0.5);

        % Significance categories
        cat_neither = valid_mask & ~x_sig & ~y_sig;
        cat_x_only = valid_mask & x_sig & ~y_sig;
        cat_y_only = valid_mask & ~x_sig & y_sig;
        cat_both = valid_mask & x_sig & y_sig;

        % Plot reference lines first (behind data)
        line([0.5 0.5], [0 1], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
        line([0 1], [0.5 0.5], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
        line([0 1], [0 1], 'Color', [0.85 0.85 0.85], 'LineStyle', '-');

        % Initialize count storage for this panel's legend
        counts = struct();
        counts.classic_neither = 0; counts.classic_x = 0; counts.classic_y = 0; counts.classic_both = 0;
        counts.intern_neither = 0; counts.intern_x = 0; counts.intern_y = 0; counts.intern_both = 0;

        % Plot each category × neuron type combination
        % Order: neither first (background), then single factors, then both (foreground)
        categories = {'neither', 'x_only', 'y_only', 'both'};
        cat_masks = {cat_neither, cat_x_only, cat_y_only, cat_both};
        cat_colors = {color_neither, color_x_only, color_y_only, color_both};

        for i_cat = 1:4
            cat_name = categories{i_cat};
            cat_mask = cat_masks{i_cat};
            cat_color = cat_colors{i_cat};

            % Classic neurons (circles 'o')
            classic_cat = cat_mask & is_classic;
            if any(classic_cat)
                scatter(x_auc(classic_cat), y_auc(classic_cat), marker_size, ...
                    'o', 'filled', 'MarkerFaceColor', cat_color, ...
                    'MarkerEdgeColor', axes_bg_color, 'LineWidth', 0.5);
            end
            counts.(['classic_' strrep(cat_name, '_only', '')]) = sum(classic_cat);

            % Interneurons (stars 'p' - pentagram)
            intern_cat = cat_mask & is_interneuron;
            if any(intern_cat)
                scatter(x_auc(intern_cat), y_auc(intern_cat), marker_size * 1.2, ...
                    'p', 'filled', 'MarkerFaceColor', cat_color, ...
                    'MarkerEdgeColor', axes_bg_color, 'LineWidth', 0.5);
            end
            counts.(['intern_' strrep(cat_name, '_only', '')]) = sum(intern_cat);
        end

        xlim([0 1]); ylim([0 1]);
        axis square;

        xlabel(sprintf('%s AUC', strrep(x_factor, '_', ' ')));
        ylabel(sprintf('%s AUC', strrep(y_factor, '_', ' ')));

        % Title with total counts
        n_classic_total = sum(valid_mask & is_classic);
        n_intern_total = sum(valid_mask & is_interneuron);
        title(sprintf('n=%d C, %d I', n_classic_total, n_intern_total), ...
            'FontWeight', 'normal');

        % Add per-panel legend with counts specific to this comparison
        x_short = capitalize_first(x_factor);
        y_short = capitalize_first(y_factor);

        % Create dummy handles for legend (in this axes)
        h_leg = gobjects(8, 1);
        leg_labels = cell(8, 1);

        % Classic entries (circles)
        h_leg(1) = scatter(NaN, NaN, legend_marker_size, 'o', 'filled', ...
            'MarkerFaceColor', color_neither, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{1} = sprintf('C: Neither (%d)', counts.classic_neither);

        h_leg(2) = scatter(NaN, NaN, legend_marker_size, 'o', 'filled', ...
            'MarkerFaceColor', color_x_only, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{2} = sprintf('C: %s (%d)', x_short, counts.classic_x);

        h_leg(3) = scatter(NaN, NaN, legend_marker_size, 'o', 'filled', ...
            'MarkerFaceColor', color_y_only, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{3} = sprintf('C: %s (%d)', y_short, counts.classic_y);

        h_leg(4) = scatter(NaN, NaN, legend_marker_size, 'o', 'filled', ...
            'MarkerFaceColor', color_both, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{4} = sprintf('C: Both (%d)', counts.classic_both);

        % Interneuron entries (stars)
        h_leg(5) = scatter(NaN, NaN, legend_marker_size * 1.2, 'p', 'filled', ...
            'MarkerFaceColor', color_neither, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{5} = sprintf('I: Neither (%d)', counts.intern_neither);

        h_leg(6) = scatter(NaN, NaN, legend_marker_size * 1.2, 'p', 'filled', ...
            'MarkerFaceColor', color_x_only, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{6} = sprintf('I: %s (%d)', x_short, counts.intern_x);

        h_leg(7) = scatter(NaN, NaN, legend_marker_size * 1.2, 'p', 'filled', ...
            'MarkerFaceColor', color_y_only, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{7} = sprintf('I: %s (%d)', y_short, counts.intern_y);

        h_leg(8) = scatter(NaN, NaN, legend_marker_size * 1.2, 'p', 'filled', ...
            'MarkerFaceColor', color_both, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{8} = sprintf('I: Both (%d)', counts.intern_both);

        legend(h_leg, leg_labels, 'Location', 'Best', 'FontSize', 6, 'Box', 'off');

    else
        % Non-SC data: use original window_roc structure with CI-based significance
        if ~isfield(aggregated_data.window_roc, epoch_name)
            continue;
        end
        if ~isfield(aggregated_data.window_roc.(epoch_name), x_factor) || ...
           ~isfield(aggregated_data.window_roc.(epoch_name), y_factor)
            continue;
        end

        x_data_all = aggregated_data.window_roc.(epoch_name).(x_factor);
        y_data_all = aggregated_data.window_roc.(epoch_name).(y_factor);

        % Concatenate across sessions
        x_auc = vertcat(x_data_all.auc);
        y_auc = vertcat(y_data_all.auc);
        x_ci = vertcat(x_data_all.ci);
        y_ci = vertcat(y_data_all.ci);

        % Valid mask
        valid_mask = ~isnan(x_auc) & ~isnan(y_auc) & ...
                     ~isnan(x_ci(:,1)) & ~isnan(y_ci(:,1));

        x_auc = x_auc(valid_mask);
        y_auc = y_auc(valid_mask);
        x_ci = x_ci(valid_mask, :);
        y_ci = y_ci(valid_mask, :);

        if isempty(x_auc)
            title(ax, sprintf('%s vs %s (No Data)', x_factor, y_factor));
            xlim([0 1]); ylim([0 1]);
            continue;
        end

        % Determine significance based on CI not containing 0.5
        x_sig = (x_ci(:,1) > 0.5) | (x_ci(:,2) < 0.5);
        y_sig = (y_ci(:,1) > 0.5) | (y_ci(:,2) < 0.5);

        % Significance categories
        cat_neither = ~x_sig & ~y_sig;
        cat_x_only = x_sig & ~y_sig;
        cat_y_only = ~x_sig & y_sig;
        cat_both = x_sig & y_sig;

        % Reference lines
        line([0.5 0.5], [0 1], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
        line([0 1], [0.5 0.5], 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
        line([0 1], [0 1], 'Color', [0.85 0.85 0.85], 'LineStyle', '-');

        % Plot each category
        if any(cat_neither)
            scatter(x_auc(cat_neither), y_auc(cat_neither), marker_size, ...
                'o', 'filled', 'MarkerFaceColor', color_neither, ...
                'MarkerEdgeColor', axes_bg_color, 'LineWidth', 0.5);
        end
        if any(cat_x_only)
            scatter(x_auc(cat_x_only), y_auc(cat_x_only), marker_size, ...
                'o', 'filled', 'MarkerFaceColor', color_x_only, ...
                'MarkerEdgeColor', axes_bg_color, 'LineWidth', 0.5);
        end
        if any(cat_y_only)
            scatter(x_auc(cat_y_only), y_auc(cat_y_only), marker_size, ...
                'o', 'filled', 'MarkerFaceColor', color_y_only, ...
                'MarkerEdgeColor', axes_bg_color, 'LineWidth', 0.5);
        end
        if any(cat_both)
            scatter(x_auc(cat_both), y_auc(cat_both), marker_size, ...
                'o', 'filled', 'MarkerFaceColor', color_both, ...
                'MarkerEdgeColor', axes_bg_color, 'LineWidth', 0.5);
        end

        xlim([0 1]); ylim([0 1]);
        axis square;

        xlabel(sprintf('%s AUC', strrep(x_factor, '_', ' ')));
        ylabel(sprintf('%s AUC', strrep(y_factor, '_', ' ')));

        % Title with counts by category
        title(sprintf('n=%d (N:%d, X:%d, Y:%d, B:%d)', length(x_auc), ...
            sum(cat_neither), sum(cat_x_only), sum(cat_y_only), sum(cat_both)), ...
            'FontWeight', 'normal', 'FontSize', 9);

        % Add per-panel legend for non-SC data
        x_short = capitalize_first(x_factor);
        y_short = capitalize_first(y_factor);

        h_leg = gobjects(4, 1);
        leg_labels = cell(4, 1);

        h_leg(1) = scatter(NaN, NaN, legend_marker_size, 'o', 'filled', ...
            'MarkerFaceColor', color_neither, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{1} = sprintf('Neither (%d)', sum(cat_neither));

        h_leg(2) = scatter(NaN, NaN, legend_marker_size, 'o', 'filled', ...
            'MarkerFaceColor', color_x_only, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{2} = sprintf('%s (%d)', x_short, sum(cat_x_only));

        h_leg(3) = scatter(NaN, NaN, legend_marker_size, 'o', 'filled', ...
            'MarkerFaceColor', color_y_only, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{3} = sprintf('%s (%d)', y_short, sum(cat_y_only));

        h_leg(4) = scatter(NaN, NaN, legend_marker_size, 'o', 'filled', ...
            'MarkerFaceColor', color_both, 'MarkerEdgeColor', axes_bg_color);
        leg_labels{4} = sprintf('Both (%d)', sum(cat_both));

        legend(h_leg, leg_labels, 'Location', 'Best', 'FontSize', 7, 'Box', 'off');
    end

    % If this is the 1st row of scatter plot panels, remove the x-tick labels
    % so they don't overlap with title information:
    if i_plot <= n_cols_scatter
        set(ax, 'XTickLabels', [])
    end
end

% General formatting
set(findall(fig2, 'Type', 'Axes'), 'TickDir', 'Out', 'LineWidth', 1, ...
    'XColor', 'k', 'YColor', 'k');

%% Save Figure 2
try
    theme('light');
catch
end

fig2_filename = fullfile(figures_dir, ...
    sprintf('window_roc_scatter_%s.pdf', brain_area_name));
pdfSave(fig2_filename, fig2.Position(3:4)/72, fig2);

end

%% ========================================================================
%  LOCAL HELPER FUNCTION
%  ========================================================================

function s = capitalize_first(str)
%CAPITALIZE_FIRST Capitalize first letter of string
    if isempty(str)
        s = str;
    else
        s = [upper(str(1)), str(2:end)];
    end
end
