%% plot_baseline_comparison.m
%
%   Generates a 2x3 summary figure from the results of the baseline
%   comparison analysis.
%
%   This script loads a session_data.mat file that is expected to contain
%   the analysis results from analyze_baseline_comparison.m.
%
% Author: Jules
% Date: 2025-09-09

%% Setup
clear; close all;
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Load Data
% This script assumes a 'session_data.mat' file is in the parent directory
% and contains the necessary analysis results.
data_file = fullfile(script_dir, '..', 'session_data.mat');
if ~exist(data_file, 'file')
    error('session_data.mat not found. Please run the analysis first.');
end
load(data_file, 'session_data');

% Extract the analysis results from the loaded data
analysis_results = session_data.analysis_results.baseline_comparison;
unique_id = session_data.metadata.unique_id;

%% Figure Setup
fig = figure('Position', [100, 100, 1200, 800]);

% Define the layout of conditions and alignments for the subplots
plot_layout = {
    {'CUE_ON', 'is_normal_dist', 'Normal Dist.'}, ...
    {'outcomeOn', 'is_norm_common', 'Common Reward'}, ...
    {'reward', 'is_common_reward_no_spe', 'Common Reward (No SPE)'}; ...
    ...
    {'CUE_ON', 'is_uniform_dist', 'Uniform Dist.'}, ...
    {'outcomeOn', 'is_norm_rare_high', 'Rare High Reward'}, ...
    {'reward', 'is_rare_high_reward_no_spe', 'Rare High Reward (No SPE)'}
    };

[n_rows, n_cols] = size(plot_layout);
axes_handles = gobjects(n_rows, n_cols);

%% Plotting Loop
for i_row = 1:n_rows
    for i_col = 1:n_cols
        % Get parameters for the current subplot
        align_name = plot_layout{i_row, i_col}{1};
        cond_name = plot_layout{i_row, i_col}{2};
        plot_title = plot_layout{i_row, i_col}{3};

        % Select the subplot
        axes_handles(i_row, i_col) = mySubPlot([n_rows, n_cols, i_col + (i_row-1)*n_cols]);
        hold on;

        % Get data for the current subplot
        sig_matrix = analysis_results.(align_name).(cond_name).sig;
        time_vector = analysis_results.(align_name).time_vector;
        n_neurons = size(sig_matrix, 1);

        % Calculate proportions
        prop_sig_up = mean(sig_matrix == 1, 1);
        prop_sig_down = -mean(sig_matrix == -1, 1); % Negative for plotting

        % Plot the proportions using barStairs
        barStairs(time_vector, prop_sig_up, 'b');
        barStairs(time_vector, prop_sig_down, 'r');

        % Formatting
        title(plot_title);
        xlim([time_vector(1), time_vector(end)]);
        ylim([-0.5, 0.5]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

        % Add labels based on position in the grid
        if i_row == n_rows
            xlabel(['Time from ' strrep(align_name, '_', ' ') ' (s)']);
        end
        if i_col == 1
            ylabel('Proportion of Neurons');
        end
    end
end

%% Final Figure Touches
% Add a legend to the first subplot
axes(axes_handles(1,1));
legend({'More Active', 'Less Active'}, 'Location', 'northwest');

% Add a main title for the whole figure
sgtitle(sprintf('Baseline vs. Post-Event Activity for: %s', unique_id), ...
    'Interpreter', 'none');

%% Check for and Plot SPE Conditions if they exist
if isfield(analysis_results.reward, 'is_common_reward_with_spe')
    fig2 = figure('Position', [100, 100, 1200, 400]);

    % Define the layout for SPE conditions
    plot_layout_spe = {
        {'reward', 'is_common_reward_with_spe', 'Common Reward (With SPE)'};
        {'reward', 'is_rare_high_reward_with_spe', 'Rare High Reward (With SPE)'}
    };

    [n_rows_spe, n_cols_spe] = size(plot_layout_spe);
    axes_handles_spe = gobjects(n_rows_spe, n_cols_spe);

    for i_col = 1:n_cols_spe % Only one row in this layout
        % Get parameters for the current subplot
        align_name = plot_layout_spe{i_col}{1};
        cond_name = plot_layout_spe{i_col}{2};
        plot_title = plot_layout_spe{i_col}{3};

        % Select the subplot
        axes_handles_spe(i_col) = mySubPlot([1, n_cols_spe, i_col]);
        hold on;

        % Get data for the current subplot
        sig_matrix = analysis_results.(align_name).(cond_name).sig;
        time_vector = analysis_results.(align_name).time_vector;
        n_neurons = size(sig_matrix, 1);

        % Calculate proportions
        prop_sig_up = mean(sig_matrix == 1, 1);
        prop_sig_down = -mean(sig_matrix == -1, 1); % Negative for plotting

        % Plot the proportions using barStairs
        barStairs(time_vector, prop_sig_up, 'b');
        barStairs(time_vector, prop_sig_down, 'r');

        % Formatting
        title(plot_title);
        xlim([time_vector(1), time_vector(end)]);
        ylim([-0.5, 0.5]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
        xlabel(['Time from ' strrep(align_name, '_', ' ') ' (s)']);

        if i_col == 1
            ylabel('Proportion of Neurons');
            legend({'More Active', 'Less Active'}, 'Location', 'northwest');
        end
    end

    sgtitle(sprintf('Baseline vs. Post-Event Activity (SPE Conditions) for: %s', unique_id), ...
        'Interpreter', 'none');
end

end
