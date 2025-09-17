%% plot_aggregated_behavior.m
%
% Generates a summary figure visualizing the proportion of sessions where
% experimental factors had a significant effect on behavioral measures.
%
% INPUTS:
%   aggregated_data   - A struct containing aggregated behavioral results.
%   brain_area_name   - A string with the name of the brain area (e.g., 'SC').
%
% Author: Jules
% Date: 2025-09-17
%
function plot_aggregated_behavior(aggregated_data, brain_area_name)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Dynamically Identify Measures and Factors
behavior_measures = fieldnames(aggregated_data.behavioral_results);
n_rows = numel(behavior_measures);

factor_names = {'Reward', 'Salience', 'Action', 'Side'};
n_cols = numel(factor_names);

if n_rows == 0
    warning('No behavioral results found in the aggregated data.');
    return;
end

%% Setup Figure
fig = figure('Position', [100, 100, 250 * n_cols, 200 * n_rows], 'Color', 'w');
h_axes = gobjects(n_rows, n_cols);

%% Plotting Loop
for i_row = 1:n_rows
    measure_name = behavior_measures{i_row};
    results_table = aggregated_data.behavioral_results.(measure_name);

    for i_col = 1:n_cols
        factor_name = factor_names{i_col};

        % --- Subplot Selection ---
        plot_idx = (i_row - 1) * n_cols + i_col;
        h_axes(i_row, i_col) = mySubPlot([n_rows, n_cols, plot_idx]);

        % --- Data Calculation ---
        % Count total unique sessions for the brain area.
        total_sessions = numel(unique(results_table.session_uuid));

        % Construct the p-value column name (e.g., 'p_Reward')
        p_value_col = ['p_', factor_name];

        % Filter for significant results
        significant_results = results_table(results_table.(p_value_col) < 0.05, :);

        % Count unique sessions that show a significant effect.
        num_significant = numel(unique(significant_results.session_uuid));

        % --- Proportion and Confidence Interval ---
        % Use binofit to get proportion and 95% CI
        [proportion, ci] = binofit(num_significant, total_sessions);

        % --- Plotting ---
        % mybarerr expects CI as [lower, upper]
        mybarerr(1, proportion, ci, 'BarWidth', 0.5);

        set(gca, 'XTick', 1, 'XTickLabel', factor_name);
        box off;
    end
end

%% Final Figure Formatting
% Use outerLims to synchronize y-axes for consistent scaling
outerLims(h_axes);

for i_row = 1:n_rows
    for i_col = 1:n_cols
        ax = h_axes(i_row, i_col);

        % Add behavioral measure as title for each row on the leftmost plot
        if i_col == 1
            title(ax, behavior_measures{i_row}, 'FontWeight', 'normal');
        end

        % Add common y-axis label to the leftmost plots
        if i_col == 1
            ylabel('Proportion of Significant Sessions');
        end

        % Remove Y-tick labels from all but the first column
        if i_col > 1
            set(ax, 'YTickLabel', []);
        end
    end
end

% Add a single, overarching title for the entire figure
sgtitle(sprintf('Aggregated Behavioral Analysis for %s', brain_area_name), ...
    'FontWeight', 'bold');

%% Save Figure
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures', 'behavior');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_behavior_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
