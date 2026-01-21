%% plot_snc_subregion_comparison.m
%
%   Generates comparison figures for SNc subregion analysis (Task 6.2).
%   This function creates two figures:
%
%   Figure 1: Grouped bar plots showing mean ROC AUC by factor and subregion,
%             organized as a grid (Rows = Epochs, Columns = Factors).
%             Error bars show standard error of the mean.
%             Statistical annotations show Wilcoxon rank-sum p-values.
%
%   Figure 2: Violin plots showing ROC AUC distributions by subregion,
%             organized as a grid (Rows = Epochs, Columns = Factors).
%             Includes individual data points overlaid on distributions.
%
% INPUTS:
%   aggregated_data - Struct containing aggregated SNc data. Must include:
%                     .neuron_metrics_table: Table with columns:
%                       - snc_subregion: 'rvmSNc' or 'cdlSNc'
%                       - is_selected: logical (putative DA neurons)
%                       - auc_{epoch}_{factor}: ROC AUC values
%
%   brain_area_name - String with brain area name (should be 'SNc').
%
% OUTPUTS:
%   None (saves PDF figures to figures/snc_subregion/ directory)
%
% Author: Claude Code
% Date: 2026-01-21
%

function plot_snc_subregion_comparison(aggregated_data, brain_area_name)

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures', 'snc_subregion');
addpath(fullfile(project_root, 'code', 'utils'));

if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

%% Input Validation
if ~strcmpi(brain_area_name, 'SNc')
    warning('plot_snc_subregion_comparison:wrong_area', ...
        'This function is designed for SNc data only. Got: %s', brain_area_name);
    return;
end

if ~isfield(aggregated_data, 'neuron_metrics_table')
    warning('plot_snc_subregion_comparison:no_table', ...
        'neuron_metrics_table is missing from the input struct.');
    return;
end

nmt = aggregated_data.neuron_metrics_table;

% Check required columns exist
required_cols = {'snc_subregion', 'is_selected'};
for i = 1:length(required_cols)
    if ~ismember(required_cols{i}, nmt.Properties.VariableNames)
        warning('plot_snc_subregion_comparison:missing_column', ...
            'Required column "%s" not found in neuron_metrics_table.', ...
            required_cols{i});
        return;
    end
end

%% Filter for Selected DA Neurons
da_neurons = nmt(nmt.is_selected, :);

% Split by subregion
rvm_mask = strcmp(da_neurons.snc_subregion, 'rvmSNc');
cdl_mask = strcmp(da_neurons.snc_subregion, 'cdlSNc');

rvm_neurons = da_neurons(rvm_mask, :);
cdl_neurons = da_neurons(cdl_mask, :);

n_rvm = height(rvm_neurons);
n_cdl = height(cdl_neurons);

if n_rvm == 0 || n_cdl == 0
    warning('plot_snc_subregion_comparison:insufficient_data', ...
        'Insufficient neurons in one or both subregions (rvmSNc: %d, cdlSNc: %d).', ...
        n_rvm, n_cdl);
    return;
end

fprintf('SNc Subregion Comparison: rvmSNc n=%d, cdlSNc n=%d\n', n_rvm, n_cdl);

%% Define Analysis Grid
epochs = {'visual', 'delay', 'perisaccade', 'postreward'};
factors = {'reward', 'probability', 'salience', 'identity'};
n_epochs = length(epochs);
n_factors = length(factors);

%% Define Colors
colors = richColors;
rvm_color = colors(6,:);   % Bright Blue for rvmSNc
cdl_color = colors(10,:);  % Bright Yellow for cdlSNc

%% ========================================================================
%  FIGURE 1: Grouped Bar Plots - Mean ROC by Subregion
%  ========================================================================

fig1 = figure('Position', [100, 100, 280 * n_factors, 200 * n_epochs], ...
    'Color', 'w', 'Name', 'SNc Subregion Comparison - Grouped Bars');
h_axes1 = gobjects(n_epochs, n_factors);

% Pre-compute statistics for all combinations
stats_data = cell(n_epochs, n_factors);

for i_epoch = 1:n_epochs
    epoch_name = epochs{i_epoch};

    for i_factor = 1:n_factors
        factor_name = factors{i_factor};
        col_name = sprintf('auc_%s_%s', epoch_name, factor_name);

        % Check if column exists
        if ~ismember(col_name, nmt.Properties.VariableNames)
            stats_data{i_epoch, i_factor} = [];
            continue;
        end

        % Extract data for each subregion
        rvm_auc = rvm_neurons.(col_name);
        cdl_auc = cdl_neurons.(col_name);

        % Remove NaN values
        rvm_auc = rvm_auc(~isnan(rvm_auc));
        cdl_auc = cdl_auc(~isnan(cdl_auc));

        % Compute statistics
        rvm_mean = mean(rvm_auc);
        cdl_mean = mean(cdl_auc);
        rvm_sem = std(rvm_auc) / sqrt(length(rvm_auc));
        cdl_sem = std(cdl_auc) / sqrt(length(cdl_auc));

        % Wilcoxon rank-sum test
        if ~isempty(rvm_auc) && ~isempty(cdl_auc)
            [p_val, ~] = ranksum(rvm_auc, cdl_auc);
        else
            p_val = NaN;
        end

        stats_data{i_epoch, i_factor} = struct(...
            'rvm_mean', rvm_mean, 'cdl_mean', cdl_mean, ...
            'rvm_sem', rvm_sem, 'cdl_sem', cdl_sem, ...
            'rvm_n', length(rvm_auc), 'cdl_n', length(cdl_auc), ...
            'p_val', p_val);
    end
end

%% Plot Figure 1
for i_epoch = 1:n_epochs
    epoch_name = epochs{i_epoch};

    for i_factor = 1:n_factors
        factor_name = factors{i_factor};

        h_axes1(i_epoch, i_factor) = mySubPlot([n_epochs, n_factors, ...
            (i_epoch-1)*n_factors + i_factor]);
        hold on;

        stats = stats_data{i_epoch, i_factor};
        if isempty(stats)
            title('No Data');
            continue;
        end

        % Bar positions
        bar_positions = [1, 2];
        bar_values = [stats.rvm_mean, stats.cdl_mean];
        bar_errors = [stats.rvm_sem, stats.cdl_sem];
        bar_colors = {rvm_color, cdl_color};

        % Plot bars
        for b = 1:2
            bar(bar_positions(b), bar_values(b), 0.6, ...
                'FaceColor', bar_colors{b}, 'EdgeColor', 'none');
        end

        % Error bars
        errorbar(bar_positions, bar_values, bar_errors, 'k', ...
            'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 6);

        % Reference line at 0.5 (chance level)
        line([0.4, 2.6], [0.5 0.5], 'Color', [0.7 0.7 0.7], ...
            'LineStyle', '--', 'LineWidth', 1);

        % Statistical annotation
        if ~isnan(stats.p_val)
            y_max = max(bar_values + bar_errors) + 0.02;
            sig_str = format_pvalue(stats.p_val);
            text(1.5, y_max + 0.02, sig_str, 'HorizontalAlignment', 'center', ...
                'FontSize', 9, 'FontWeight', 'normal');

            % Significance bracket
            if stats.p_val < 0.05
                line([1, 2], [y_max, y_max], 'Color', 'k', 'LineWidth', 1);
            end
        end

        % Axis formatting
        xlim([0.4, 2.6]);
        ylim([0.35, 0.75]);
        set(gca, 'XTick', [1, 2], 'XTickLabel', {'rvmSNc', 'cdlSNc'});
    end
end

%% Figure 1 Formatting
for i_epoch = 1:n_epochs
    for i_factor = 1:n_factors
        h_ax = h_axes1(i_epoch, i_factor);
        if ~isgraphics(h_ax); continue; end

        % Column headers (top row only)
        if i_epoch == 1
            title(h_ax, strrep(factors{i_factor}, '_', ' '), ...
                'FontWeight', 'bold');
        end

        % Row labels (first column only)
        if i_factor == 1
            ylabel(h_ax, sprintf('%s\nAUC', strrep(epochs{i_epoch}, '_', ' ')));
        else
            set(h_ax, 'YTickLabel', []);
        end

        % Remove x-axis labels except bottom row
        if i_epoch < n_epochs
            set(h_ax, 'XTickLabel', []);
        end
    end
end

% General formatting
all_valid_axes = h_axes1(isgraphics(h_axes1));
set(all_valid_axes, 'TickDir', 'Out', 'Color', 'none', ...
    'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, 'FontSize', 9);

% Add legend
legend_ax = h_axes1(1, n_factors);
axes(legend_ax);
h_leg = legend({'rvmSNc', 'cdlSNc'}, 'Location', 'northeast', ...
    'FontSize', 8, 'Box', 'off');

%% Save Figure 1
fig1_filename = fullfile(figures_dir, 'snc_subregion_grouped_bars.pdf');
pdfSave(fig1_filename, fig1.Position(3:4)/72, fig1);

%% ========================================================================
%  FIGURE 2: Violin Plots - ROC Distributions by Subregion
%  ========================================================================

fig2 = figure('Position', [100, 100, 280 * n_factors, 220 * n_epochs], ...
    'Color', 'w', 'Name', 'SNc Subregion Comparison - Violin Plots');
h_axes2 = gobjects(n_epochs, n_factors);

for i_epoch = 1:n_epochs
    epoch_name = epochs{i_epoch};

    for i_factor = 1:n_factors
        factor_name = factors{i_factor};
        col_name = sprintf('auc_%s_%s', epoch_name, factor_name);

        h_axes2(i_epoch, i_factor) = mySubPlot([n_epochs, n_factors, ...
            (i_epoch-1)*n_factors + i_factor]);
        hold on;

        % Check if column exists
        if ~ismember(col_name, nmt.Properties.VariableNames)
            title('No Data');
            continue;
        end

        % Extract data
        rvm_auc = rvm_neurons.(col_name);
        cdl_auc = cdl_neurons.(col_name);
        rvm_auc = rvm_auc(~isnan(rvm_auc));
        cdl_auc = cdl_auc(~isnan(cdl_auc));

        % Plot violins
        plot_violin(1, rvm_auc, rvm_color, 0.35);
        plot_violin(2, cdl_auc, cdl_color, 0.35);

        % Add individual data points with jitter
        jitter_width = 0.1;
        scatter(1 + jitter_width * (rand(length(rvm_auc), 1) - 0.5), ...
            rvm_auc, 12, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
        scatter(2 + jitter_width * (rand(length(cdl_auc), 1) - 0.5), ...
            cdl_auc, 12, 'k', 'filled', 'MarkerFaceAlpha', 0.3);

        % Add median lines
        line([0.7, 1.3], [median(rvm_auc), median(rvm_auc)], ...
            'Color', 'k', 'LineWidth', 2);
        line([1.7, 2.3], [median(cdl_auc), median(cdl_auc)], ...
            'Color', 'k', 'LineWidth', 2);

        % Reference line at 0.5 (chance level)
        line([0.4, 2.6], [0.5 0.5], 'Color', [0.7 0.7 0.7], ...
            'LineStyle', '--', 'LineWidth', 1);

        % Statistical annotation
        stats = stats_data{i_epoch, i_factor};
        if ~isempty(stats) && ~isnan(stats.p_val)
            y_max = max([rvm_auc; cdl_auc]) + 0.05;
            sig_str = format_pvalue(stats.p_val);
            text(1.5, y_max + 0.03, sig_str, 'HorizontalAlignment', 'center', ...
                'FontSize', 9, 'FontWeight', 'normal');

            % Significance bracket
            if stats.p_val < 0.05
                line([1, 2], [y_max, y_max], 'Color', 'k', 'LineWidth', 1);
            end
        end

        % Axis formatting
        xlim([0.4, 2.6]);
        ylim([0.2, 1.0]);
        set(gca, 'XTick', [1, 2], 'XTickLabel', {'rvmSNc', 'cdlSNc'});
    end
end

%% Figure 2 Formatting
for i_epoch = 1:n_epochs
    for i_factor = 1:n_factors
        h_ax = h_axes2(i_epoch, i_factor);
        if ~isgraphics(h_ax); continue; end

        % Column headers (top row only)
        if i_epoch == 1
            title(h_ax, strrep(factors{i_factor}, '_', ' '), ...
                'FontWeight', 'bold');
        end

        % Row labels (first column only)
        if i_factor == 1
            ylabel(h_ax, sprintf('%s\nAUC', strrep(epochs{i_epoch}, '_', ' ')));
        else
            set(h_ax, 'YTickLabel', []);
        end

        % Remove x-axis labels except bottom row
        if i_epoch < n_epochs
            set(h_ax, 'XTickLabel', []);
        end
    end
end

% General formatting
all_valid_axes = h_axes2(isgraphics(h_axes2));
set(all_valid_axes, 'TickDir', 'Out', 'Color', 'none', ...
    'XColor', 'k', 'YColor', 'k', 'LineWidth', 1, 'FontSize', 9);

%% Save Figure 2
fig2_filename = fullfile(figures_dir, 'snc_subregion_violin_plots.pdf');
pdfSave(fig2_filename, fig2.Position(3:4)/72, fig2);

%% Print Summary Statistics Table
fprintf('\n=== SNc Subregion Comparison Summary ===\n');
fprintf('%-12s %-12s | rvmSNc (n=%d) | cdlSNc (n=%d) | p-value\n', ...
    'Epoch', 'Factor', n_rvm, n_cdl);
fprintf('%s\n', repmat('-', 1, 70));

for i_epoch = 1:n_epochs
    for i_factor = 1:n_factors
        stats = stats_data{i_epoch, i_factor};
        if isempty(stats); continue; end

        fprintf('%-12s %-12s | %.3f +/- %.3f | %.3f +/- %.3f | %s\n', ...
            epochs{i_epoch}, factors{i_factor}, ...
            stats.rvm_mean, stats.rvm_sem, ...
            stats.cdl_mean, stats.cdl_sem, ...
            format_pvalue(stats.p_val));
    end
end

fprintf('\nFigures saved to: %s\n', figures_dir);

end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function sig_str = format_pvalue(p)
%FORMAT_PVALUE Format p-value for display with significance markers
    if isnan(p)
        sig_str = 'n/a';
    elseif p < 0.001
        sig_str = 'p < 0.001 ***';
    elseif p < 0.01
        sig_str = sprintf('p = %.3f **', p);
    elseif p < 0.05
        sig_str = sprintf('p = %.3f *', p);
    else
        sig_str = sprintf('p = %.2f', p);
    end
end

function plot_violin(x_pos, data, color, width)
%PLOT_VIOLIN Simple violin plot using kernel density estimation
%   x_pos - center x position for the violin
%   data - vector of values
%   color - RGB triplet for fill color
%   width - half-width of the violin at maximum density

    if isempty(data) || length(data) < 3
        return;
    end

    % Compute kernel density estimate
    [f, xi] = ksdensity(data, 'NumPoints', 100);

    % Normalize density to fit within width
    f_norm = f / max(f) * width;

    % Create violin shape (symmetric around x_pos)
    x_violin = [x_pos + f_norm, x_pos - flip(f_norm)];
    y_violin = [xi, flip(xi)];

    % Plot filled violin
    fill(x_violin, y_violin, color, 'EdgeColor', color * 0.7, ...
        'FaceAlpha', 0.5, 'LineWidth', 1);
end
