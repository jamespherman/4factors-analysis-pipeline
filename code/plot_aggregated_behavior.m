%% plot_aggregated_behavior.m
%
% Generates a summary figure visualizing the proportion of sessions
% with significant behavioral effects, based on the two-model
% (Image vs. Bullseye trials) analysis. The figure is structured as a
% 3x2 grid: Rows for behavioral measures and columns for main vs.
% interaction effects.
%
% INPUTS:
%   aggregated_data   - A struct with aggregated behavioral results.
%   brain_area_name   - A string with the name of the brain area.
%
% Author: Jules
% Date: 2025-09-21
%
function plot_aggregated_behavior(aggregated_data, brain_area_name)

%% Setup
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

measures = {'reaction_time', 'peak_velocity', 'endpoint_error'};
model_types = {'imagetrials', 'bullseyetrials'};
n_rows = numel(measures);
n_cols = 2; % Main Effects, Interaction Effects

fig = figure('Position', [100, 100, 800, 900], 'Color', 'w');
h_axes = gobjects(n_rows, n_cols);

colors = richColors(2); % Colors for Image and Bullseye models

%% Plotting Loop
for i_row = 1:n_rows
    measure = measures{i_row};

    % --- Column 1: Main Effects ---
    ax_main = mySubPlot([n_rows, n_cols, (i_row - 1) * n_cols + 1]);
    h_axes(i_row, 1) = ax_main;
    hold on;

    main_effects_img = {'reward', 'probability', 'identity'};
    main_effects_bull = {'reward', 'probability', 'saliency'};
    labels_main = {'Reward', 'Probability', 'Identity/Saliency'};

    proportions = nan(2, 3); % 2 models, 3 effects

    % Image Trials
    table_img = aggregated_data.behavioral_results.([measure '_imagetrials']);
    sessions_img = unique(table_img.session_id);
    for i = 1:3
        sig_sessions = unique(table_img.session_id( ...
            strcmp(table_img.Effect, main_effects_img{i}) & table_img.pValue < 0.05));
        proportions(1, i) = numel(sig_sessions) / numel(sessions_img);
    end

    % Bullseye Trials
    table_bull = aggregated_data.behavioral_results.([measure '_bullseyetrials']);
    sessions_bull = unique(table_bull.session_id);
    for i = 1:3
        sig_sessions = unique(table_bull.session_id( ...
            strcmp(table_bull.Effect, main_effects_bull{i}) & table_bull.pValue < 0.05));
        proportions(2, i) = numel(sig_sessions) / numel(sessions_bull);
    end

    b_main = bar(proportions', 'grouped');
    b_main(1).FaceColor = colors(1, :);
    b_main(2).FaceColor = colors(2, :);
    set(ax_main, 'XTickLabel', labels_main);
    title(ax_main, 'Main Effects');

    % --- Column 2: Interaction Effects ---
    ax_int = mySubPlot([n_rows, n_cols, (i_row - 1) * n_cols + 2]);
    h_axes(i_row, 2) = ax_int;
    hold on;

    int_effects_img = {'reward:probability', 'reward:identity', 'probability:identity'};
    int_effects_bull = {'reward:probability', 'reward:saliency', 'probability:saliency'};
    labels_int = {'R x P', 'R x I/S', 'P x I/S'};

    proportions_int = nan(2, 3); % 2 models, 3 effects

    % Image Trials
    for i = 1:3
        sig_sessions = unique(table_img.session_id( ...
            strcmp(table_img.Effect, int_effects_img{i}) & table_img.pValue < 0.05));
        proportions_int(1, i) = numel(sig_sessions) / numel(sessions_img);
    end

    % Bullseye Trials
    for i = 1:3
        sig_sessions = unique(table_bull.session_id( ...
            strcmp(table_bull.Effect, int_effects_bull{i}) & table_bull.pValue < 0.05));
        proportions_int(2, i) = numel(sig_sessions) / numel(sessions_bull);
    end

    b_int = bar(proportions_int', 'grouped');
    b_int(1).FaceColor = colors(1, :);
    b_int(2).FaceColor = colors(2, :);
    set(ax_int, 'XTickLabel', labels_int);
    title(ax_int, 'Interaction Effects');
end

%% Final Figure Formatting
linkaxes(h_axes(:), 'y');
ylim(h_axes(1,1), [0, 1]);

for i_row = 1:n_rows
    ylabel(h_axes(i_row, 1), measures{i_row}, 'FontWeight', 'bold');
    set(h_axes(i_row, 2), 'YTickLabel', []);
end

legend(h_axes(1,2), {'Image Trials', 'Bullseye Trials'}, 'Location', 'northeast');
sgtitle(sprintf('Aggregated Behavioral Analysis for %s', brain_area_name), ...
    'FontWeight', 'bold');

%% Save Figure
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures', 'behavior');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_behavior_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
