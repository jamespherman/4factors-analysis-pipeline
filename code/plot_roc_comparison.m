%% plot_roc_comparison.m
%
%   Generates a 1x3 summary figure from the results of the between-condition
%   ROC comparison analysis.
%
%   This script loads a session_data.mat file that is expected to contain
%   the analysis results from analyze_roc_comparison.m.
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
analysis_results = session_data.analysis_results.roc_comparison;
unique_id = session_data.metadata.unique_id;

%% Figure Setup
fig = figure('Position', [100, 100, 1500, 500]);

% The results struct contains the names of the comparisons
comparison_names = fieldnames(analysis_results);
n_comparisons = numel(comparison_names);

%% Plotting Loop
for i_comp = 1:n_comparisons
    comp_name = comparison_names{i_comp};
    comp_data = analysis_results.(comp_name);

    % Select the subplot
    mySubPlot([1, n_comparisons, i_comp]);
    hold on;

    % Get data for the current subplot
    sig_matrix = comp_data.sig;
    time_vector = comp_data.time_vector;
    cond_names = comp_data.cond_names;

    % Clean up condition names for the legend
    legend_name_1 = strrep(cond_names{1}, 'is_', '');
    legend_name_1 = strrep(legend_name_1, '_', ' ');
    legend_name_2 = strrep(cond_names{2}, 'is_', '');
    legend_name_2 = strrep(legend_name_2, '_', ' ');

    % Calculate proportions of neurons preferring each condition
    % arrayROC(cond1, cond2): sig=1 means cond2>cond1, sig=-1 means cond1>cond2
    prop_cond2_pref = mean(sig_matrix == 1, 1);
    prop_cond1_pref = -mean(sig_matrix == -1, 1); % Negative for plotting

    % Plot the proportions using barStairs
    h1 = barStairs(time_vector, prop_cond2_pref, 'g');
    h2 = barStairs(time_vector, prop_cond1_pref, 'm');

    % Formatting
    % Extract alignment event from comparison name for the title
    title_text = strrep(comp_name, '_', ' ');
    title(title_text);
    xlim([time_vector(1), time_vector(end)]);
    ylim([-0.5, 0.5]);
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
    xlabel('Time (s)');

    % Add legend to each subplot
    legend([h1, h2], {legend_name_2, legend_name_1}, 'Location', 'northwest');

    % Add Y-axis label to the first plot only
    if i_comp == 1
        ylabel('Proportion of Neurons');
    end
end

%% Final Figure Touches
sgtitle(sprintf('Between-Condition ROC Comparison for: %s', unique_id), ...
    'Interpreter', 'none');

end
