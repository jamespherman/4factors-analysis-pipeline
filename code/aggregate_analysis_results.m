%% aggregate_analysis_results.m
%
% Systematically gathers all pre-computed single-session analysis results
% and pools them by brain area. This script now relies on a predefined
% analysis plan, leveraging the standardized data structure.
%
% OUTPUT:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-13 (Refactored)
%

function [aggregated_sc_data, aggregated_snc_data] = aggregate_analysis_results

% Start timer and define in-line function to give user feedback:
tic;
giveFeed = @(x)disp([num2str(round(toc, 1)) 's - ' x]);

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
addpath(project_root);

%% Load Manifest and Analysis Plan
giveFeed('Loading session manifest and analysis plan...');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);
[~, ~, analysis_plan] = define_tokens_task_conditions(); % Get the canonical plan
baseline_events = {'CUE_ON', 'outcomeOn', 'reward'}; % This is hard-coded in the analysis fn
giveFeed('Manifest and plan loaded.');

%% Initialize Aggregated Data Structures
brain_areas = {'SC', 'SNc'};
aggregated_data_by_area = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i_area = 1:length(brain_areas)
    area = brain_areas{i_area};
    agg_data = struct();
    agg_data.session_id = {};

    % Initialize based on the known analysis plan
    % 1. ROC Comparisons
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);
        agg_data.roc_comparison.(comp.event).(comp.name).sig = [];
    end
    % 2. Baseline Comparisons
    for j = 1:length(analysis_plan.baseline_plan)
        comp = analysis_plan.baseline_plan(j);
        for i_event = 1:length(baseline_events)
            event_name = baseline_events{i_event};
            agg_data.baseline_comparison.(event_name).(comp.name).sig = [];
        end
    end
    % 3. ANOVA Results
    if analysis_plan.anova_plan.run
        for i_event = 1:length(analysis_plan.anova_plan.events)
            event_name = analysis_plan.anova_plan.events{i_event};
            for i_p = 1:length(analysis_plan.anova_plan.p_value_names)
                p_name = analysis_plan.anova_plan.p_value_names{i_p};
                agg_data.anova_results.(event_name).(p_name) = [];
            end
        end
    end
    aggregated_data_by_area(area) = agg_data;
end

%% Aggregation Loop
complete_sessions = manifest(strcmp(manifest.analysis_status, 'complete'), :);
nSessions = height(complete_sessions);
giveFeed(sprintf('Starting aggregation for %d completed sessions.', nSessions));

for i_session = 1:nSessions
    session_id = complete_sessions.unique_id{i_session};
    brain_area = complete_sessions.brain_area{i_session};
    giveFeed(sprintf('Processing session %d/%d: %s', i_session, nSessions, session_id));

    session_data_path = fullfile(findOneDrive(), 'Neuronal Data Analysis', ...
        session_id, [session_id '_session_data.mat']);

    if ~exist(session_data_path, 'file')
        warning('aggregate_analysis_results:fileNotFound', ...
                'Could not find session_data.mat for %s. Skipping.', session_id);
        continue;
    end

    data = load(session_data_path, 'session_data');
    session_data = data.session_data;

    if ~isfield(session_data, 'analysis')
        warning('aggregate_analysis_results:noAnalysis', ...
            'No analysis field in data for %s. Skipping.', session_id);
        continue;
    end

    n_neurons = size(session_data.analysis.selected_neurons, 1);

    % Get the correct aggregated struct for the current brain area
    agg_data = aggregated_data_by_area(brain_area);

    % Append session ID
    agg_data.session_id = [agg_data.session_id; repmat({session_id}, n_neurons, 1)];

    % --- 1. Aggregate ROC Comparison Results ---
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);
        event = comp.event;
        name = comp.name;

        data_to_append = nan(n_neurons, 1); % Default to NaN
        if isfield(session_data.analysis, 'roc_comparison') && ...
           isfield(session_data.analysis.roc_comparison, event) && ...
           isfield(session_data.analysis.roc_comparison.(event), name) && ...
           isfield(session_data.analysis.roc_comparison.(event).(name), 'sig')

            data_to_append = session_data.analysis.roc_comparison.(event).(name).sig;
            if ~isempty(data_to_append) && ...
               ~isfield(agg_data.roc_comparison.(event).(name), 'time_vector')
                agg_data.roc_comparison.(event).(name).time_vector = ...
                    session_data.analysis.roc_comparison.(event).(name).time_vector;
            end
        end
        agg_data.roc_comparison.(event).(name).sig = ...
            [agg_data.roc_comparison.(event).(name).sig; data_to_append];
    end

    % --- 2. Aggregate Baseline Comparison Results ---
    for j = 1:length(analysis_plan.baseline_plan)
        comp = analysis_plan.baseline_plan(j);
        name = comp.name;
        for i_event = 1:length(baseline_events)
            event = baseline_events{i_event};

            data_to_append = nan(n_neurons, 1); % Default to NaN
            if isfield(session_data.analysis, 'baseline_comparison') && ...
               isfield(session_data.analysis.baseline_comparison, event) && ...
               isfield(session_data.analysis.baseline_comparison.(event), name) && ...
               isfield(session_data.analysis.baseline_comparison.(event).(name), 'sig')

                data_to_append = session_data.analysis.baseline_comparison.(event).(name).sig;
                if ~isempty(data_to_append) && ...
                   ~isfield(agg_data.baseline_comparison.(event).(name), 'time_vector')
                    agg_data.baseline_comparison.(event).(name).time_vector = ...
                        session_data.analysis.baseline_comparison.(event).(name).time_vector;
                end
            end
            agg_data.baseline_comparison.(event).(name).sig = ...
                [agg_data.baseline_comparison.(event).(name).sig; data_to_append];
        end
    end

    % --- 3. Aggregate ANOVA Results ---
    if analysis_plan.anova_plan.run && isfield(session_data.analysis, 'anova_results')
        for i_event = 1:length(analysis_plan.anova_plan.events)
            event = analysis_plan.anova_plan.events{i_event};
            for i_p = 1:length(analysis_plan.anova_plan.p_value_names)
                p_name = analysis_plan.anova_plan.p_value_names{i_p};

                data_to_append = nan(n_neurons, 1); % Default to NaN
                if isfield(session_data.analysis.anova_results, event) && ...
                   isfield(session_data.analysis.anova_results.(event), p_name)
                    data_to_append = session_data.analysis.anova_results.(event).(p_name);
                end
                agg_data.anova_results.(event).(p_name) = ...
                    [agg_data.anova_results.(event).(p_name); data_to_append];
            end
        end
    end

    % Update the map with the modified struct
    aggregated_data_by_area(brain_area) = agg_data;
end

%% Finalize and Save
giveFeed('Saving aggregated data to file...');
aggregated_sc_data = aggregated_data_by_area('SC');
aggregated_snc_data = aggregated_data_by_area('SNc');

saveFileName = fullfile(project_root, 'data', 'processed', ...
    'aggregated_analysis_data.mat');
save(saveFileName, 'aggregated_sc_data', 'aggregated_snc_data', '-v7.3');
giveFeed('Aggregation complete.');

end
