%% plot_aggregated_behavior.m
%
%   Generates a summary figure visualizing the proportion of sessions with
%   significant behavioral effects. The output adheres to the specifications
%   in `docs/plotting_requirements.md`.
%
%   The figure is structured with behavioral measures as rows and effect types
%   (Main vs. Interaction) as columns.
%
% INPUTS:
%
%   aggregated_data - A struct containing aggregated data for one brain area.
%                     This script requires the `behavioral_results` field,
%                     which is a struct of tables, one for each analysis name
%                     (e.g., 'reaction_time_imagetrials'). Each table must
%                     be in a "tidy" format and contain 'Effect', 'pValue',
%                     and 'session_id' columns.
%
%   brain_area_name - A string with the name of the brain area (e.g., 'SC').
%   analysis_plan   - The analysis plan struct, which is the source of truth
%                     for behavioral measures and model definitions.
%
% Author: Jules
% Date: 2025-09-24
%
function plot_aggregated_behavior(aggregated_data, brain_area_name, analysis_plan)

%% Setup
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

if ~isfield(aggregated_data, 'behavioral_results')
    warning('plot_aggregated_behavior:no_data', 'No behavioral data found.');
    return;
end

%% Dynamically Determine Plot Layout from Analysis Plan
behavior_plan = analysis_plan.behavior_plan;
measures = {behavior_plan.name};
models = analysis_plan.behavior_models;
model_labels = {models.label};

% Find all unique main and interaction effects across all models
all_main_effects = {};
all_interaction_effects = {};
for i = 1:length(models)
    all_main_effects = union(all_main_effects, models(i).main_effects);
    all_interaction_effects = union(all_interaction_effects, models(i).interaction_effects);
end

n_rows = numel(measures);
n_cols = 2; % Main Effects, Interaction Effects

fig = figure('Position', [100, 100, 400 * n_cols, 300 * n_rows], 'Color', 'w');
h_axes = gobjects(n_rows, n_cols);
colors = richColors(length(models));

%% Plotting Loop
for i_row = 1:n_rows
    measure = measures{i_row};

    % --- Column 1: Main Effects ---
    ax_main = mySubPlot([n_rows, n_cols, (i_row - 1) * n_cols + 1]);
    h_axes(i_row, 1) = ax_main;
    hold on;

    proportions_main = nan(length(models), length(all_main_effects));
    for i_model = 1:length(models)
        model_def = models(i_model);
        analysis_name = [measure, '_', model_def.name];
        if ~isfield(aggregated_data.behavioral_results, analysis_name)
            continue;
        end

        tbl = aggregated_data.behavioral_results.(analysis_name);
        all_sessions = unique(tbl.session_id);
        n_total_sessions = numel(all_sessions);
        if n_total_sessions == 0, continue; end

        for i_effect = 1:length(all_main_effects)
            effect_name = all_main_effects{i_effect};
            if ismember(effect_name, model_def.main_effects)
                sig_sessions = unique(tbl.session_id(strcmp(tbl.Effect, effect_name) & tbl.pValue < 0.05));
                proportions_main(i_model, i_effect) = numel(sig_sessions) / n_total_sessions;
            end
        end
    end

    b_main = bar(proportions_main', 'grouped');
    for i = 1:length(b_main), b_main(i).FaceColor = colors(i,:); end
    set(ax_main, 'XTickLabel', all_main_effects, 'XTickLabelRotation', 30);
    title(ax_main, 'Main Effects');

    % --- Column 2: Interaction Effects ---
    ax_int = mySubPlot([n_rows, n_cols, (i_row - 1) * n_cols + 2]);
    h_axes(i_row, 2) = ax_int;
    hold on;

    proportions_int = nan(length(models), length(all_interaction_effects));
    for i_model = 1:length(models)
        model_def = models(i_model);
        analysis_name = [measure, '_', model_def.name];
        if ~isfield(aggregated_data.behavioral_results, analysis_name)
            continue;
        end

        tbl = aggregated_data.behavioral_results.(analysis_name);
        all_sessions = unique(tbl.session_id);
        n_total_sessions = numel(all_sessions);
        if n_total_sessions == 0, continue; end

        for i_effect = 1:length(all_interaction_effects)
            effect_name = all_interaction_effects{i_effect};
             if ismember(effect_name, model_def.interaction_effects)
                sig_sessions = unique(tbl.session_id(strcmp(tbl.Effect, effect_name) & tbl.pValue < 0.05));
                proportions_int(i_model, i_effect) = numel(sig_sessions) / n_total_sessions;
            end
        end
    end

    b_int = bar(proportions_int', 'grouped');
    for i = 1:length(b_int), b_int(i).FaceColor = colors(i,:); end
    set(ax_int, 'XTickLabel', strrep(all_interaction_effects, ':', ' x '), 'XTickLabelRotation', 30);
    title(ax_int, 'Interaction Effects');
end

%% Final Figure Formatting
linkaxes(h_axes(isgraphics(h_axes)), 'y');
ylim(h_axes(1,1), [0, 1]); % Set Y-axis from 0 to 1

for i_row = 1:n_rows
    if isgraphics(h_axes(i_row, 1))
        ylabel(h_axes(i_row, 1), strrep(measures{i_row}, '_', ' '), 'FontWeight', 'bold');
    end
    if isgraphics(h_axes(i_row, 2))
        set(h_axes(i_row, 2), 'YTickLabel', []);
    end
end

legend(h_axes(1,2), model_labels, 'Location', 'northeast', 'Interpreter', 'none');
sgtitle(sprintf('Aggregated Behavioral Analysis for %s', brain_area_name), ...
    'FontWeight', 'bold', 'Interpreter', 'none');

%% Save Figure
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_behavioral_effects_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
