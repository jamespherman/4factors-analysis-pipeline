%% plot_aggregated_behavior.m
%
% Generates a summary figure visualizing the proportion of sessions
% with significant behavioral effects. The script is plan-driven,
% deriving its layout and data from the analysis_plan and the aggregated
% results file.
%
% INPUTS:
%   aggregated_data   - A struct with aggregated behavioral results.
%   brain_area_name   - A string with the name of the brain area.
%   analysis_plan     - The analysis plan struct from define_task_conditions.
%
% Author: Jules
% Date: 2025-09-25 (Refactored for dynamic, plan-driven layout)
%
function plot_aggregated_behavior(aggregated_data, brain_area_name, analysis_plan)

%% Setup
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Dynamically Determine Plot Layout from Analysis Plan
behavior_plan = analysis_plan.behavior_plan;

% 1. Get a unique list of behavioral measures (e.g., 'reaction_time')
all_analysis_names = {behavior_plan.name};
measure_names = unique(cellfun(@(x) extractBefore(x, '_'), ...
    all_analysis_names, 'UniformOutput', false));
n_rows = numel(measure_names);

% 2. Get a unique list of model names (e.g., 'imagetrials') and their labels
model_names = unique(cellfun(@(x) extractAfter(x, '_'), ...
    all_analysis_names, 'UniformOutput', false));
n_models = numel(model_names);
model_labels = cell(1, n_models);
for i = 1:n_models
    if contains(model_names{i}, 'image')
        model_labels{i} = 'Image Trials';
    elseif contains(model_names{i}, 'bullseye')
        model_labels{i} = 'Bullseye Trials';
    else
        model_labels{i} = model_names{i};
    end
end

% 3. Find all unique main and interaction effects from the aggregated data
try
all_effects = {};
for i = 1:length(all_analysis_names)
    if isfield(aggregated_data.behavioral_results, all_analysis_names{i})
        tbl = aggregated_data.behavioral_results.(all_analysis_names{i});
        all_effects = union(all_effects, tbl.Effect);
    end
end
catch me
    keyboard
end
all_effects = setdiff(all_effects, '(Intercept)'); % Remove intercept
main_effects = all_effects(~contains(all_effects, ':'));
interaction_effects = all_effects(contains(all_effects, ':'));

% --- Figure and Plotting Setup ---
n_cols = 2; % Main Effects, Interaction Effects
fig = figure('Position', [100, 100, 800, 300 * n_rows], 'Color', 'w');
h_axes = gobjects(n_rows, n_cols);
colors = richColors(n_models); % Use the project's standard palette

%% Plotting Loop
for i_row = 1:n_rows
    measure = measure_names{i_row};

    % --- Column 1: Main Effects ---
    ax_main = mySubPlot([n_rows, n_cols, (i_row - 1) * n_cols + 1]);
    h_axes(i_row, 1) = ax_main;
    hold on;
    
    proportions_main = nan(n_models, length(main_effects));

    for i_model = 1:n_models
        model_name = model_names{i_model};
        analysis_name = [measure, '_', model_name];

        if ~isfield(aggregated_data.behavioral_results, analysis_name)
            continue;
        end
        
        tbl = aggregated_data.behavioral_results.(analysis_name);
        all_sessions = unique(tbl.session_id);
        n_total_sessions = numel(all_sessions);

        if n_total_sessions == 0, continue; end

        for i_effect = 1:length(main_effects)
            effect_name = main_effects{i_effect};
            
            plan_idx = find(strcmp({behavior_plan.name}, analysis_name));
            if ismember(effect_name, behavior_plan(plan_idx).factors)
                sig_sessions = unique(tbl.session_id(strcmp(tbl.Effect, effect_name) & tbl.pValue < 0.05));
                proportions_main(i_model, i_effect) = numel(sig_sessions) / n_total_sessions;
            end
        end
    end

    b_main = bar(proportions_main', 'grouped');
    for i = 1:length(b_main), b_main(i).FaceColor = colors(i,:); end
    set(ax_main, 'XTickLabel', main_effects, 'XTickLabelRotation', 30);
    title(ax_main, 'Main Effects');

    % --- Column 2: Interaction Effects ---
    ax_int = mySubPlot([n_rows, n_cols, (i_row - 1) * n_cols + 2]);
    h_axes(i_row, 2) = ax_int;
    hold on;

    proportions_int = nan(n_models, length(interaction_effects));

    for i_model = 1:n_models
        model_name = model_names{i_model};
        analysis_name = [measure, '_', model_name];

        if ~isfield(aggregated_data.behavioral_results, analysis_name)
            continue;
        end

        tbl = aggregated_data.behavioral_results.(analysis_name);
        all_sessions = unique(tbl.session_id);
        n_total_sessions = numel(all_sessions);

        if n_total_sessions == 0, continue; end

        for i_effect = 1:length(interaction_effects)
            effect_name = interaction_effects{i_effect};
            
            plan_idx = find(strcmp({behavior_plan.name}, analysis_name));
            model_factors = behavior_plan(plan_idx).factors;
            interaction_factors = strsplit(effect_name, ':');
            
            if all(ismember(interaction_factors, model_factors))
                sig_sessions = unique(tbl.session_id(strcmp(tbl.Effect, effect_name) & tbl.pValue < 0.05));
                proportions_int(i_model, i_effect) = numel(sig_sessions) / n_total_sessions;
            end
        end
    end
    
    b_int = bar(proportions_int', 'grouped');
    for i = 1:length(b_int), b_int(i).FaceColor = colors(i,:); end
    set(ax_int, 'XTickLabel', strrep(interaction_effects, ':', ' x '), 'XTickLabelRotation', 30);
    title(ax_int, 'Interaction Effects');
end

%% Final Figure Formatting
linkaxes(h_axes(isgraphics(h_axes)), 'y');
ylim(h_axes(1,1), [0, 1]);

for i_row = 1:n_rows
    ylabel(h_axes(i_row, 1), strrep(measure_names{i_row}, '_', ' '), 'FontWeight', 'bold');
    set(h_axes(i_row, 2), 'YTickLabel', []);
end

legend(h_axes(1,2), model_labels, 'Location', 'northeast');
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