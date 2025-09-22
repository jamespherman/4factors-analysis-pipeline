%% aggregate_analysis_results.m
%
% Gathers all single-session analysis results and pools them by
% brain area. This script is aligned with the two-model analysis
% plan defined in `define_task_conditions.m`, handling the nested
% output structures for Image and Bullseye trial analyses.
%
% OUTPUT:
%   aggregated_sc_data  - A struct with aggregated data for SC.
%   aggregated_snc_data - A struct with aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-21 (Refactored for two-model approach)
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
[~, analysis_plan] = define_task_conditions(); % Get the canonical plan
baseline_events = {'CUE_ON', 'outcomeOn', 'reward'}; % This is hard-coded in the analysis fn
giveFeed('Manifest and plan loaded.');

%% Initialize Aggregated Data Structures
% This section initializes the master data structures based on the canonical
% analysis_plan. This replaces any old logic that dynamically discovered
% the analysis structure.
brain_areas = {'SC', 'SNc'};
aggregated_data_by_area = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i_area = 1:length(brain_areas)
    area = brain_areas{i_area};
    agg_data = struct();

    % This field will be populated with neuron-level session IDs
    agg_data.session_id = {};

    % 1. Initialize ROC Comparison Results
    agg_data.roc_comparison = struct();
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);
        agg_data.roc_comparison.(comp.event).(comp.name).sig = [];
    end

    % 2. Initialize Baseline Comparison Results
    agg_data.baseline_comparison = struct();
    for j = 1:length(analysis_plan.baseline_plan)
        comp = analysis_plan.baseline_plan(j);
        for i_event = 1:length(baseline_events)
            event_name = baseline_events{i_event};
            agg_data.baseline_comparison.(event_name).(comp.name).sig = [];
        end
    end

    % 3. Initialize ANOVA Results (for Two-Model Structure)
    agg_data.anova_results = struct();
    for j = 1:length(analysis_plan.anova_plan)
        plan_item = analysis_plan.anova_plan(j);
        if plan_item.run
            % The results from analyze_anova are not nested by event
            % but are stored directly under the analysis name.
            agg_data.anova_results.(plan_item.name) = struct();
        end
    end

    % 4. Initialize Behavioral Results (for Two-Model Structure)
    agg_data.behavioral_results = struct();
    for j = 1:length(analysis_plan.behavior_plan)
        plan_item = analysis_plan.behavior_plan(j);
        agg_data.behavioral_results.(plan_item.name) = table();
    end

    % 5. Initialize Population Decoding Results
    agg_data.population_decoding = struct();
    for j = 1:length(analysis_plan.decoding_plan.testing_plan)
        test_name = analysis_plan.decoding_plan.testing_plan(j).test_name;
        agg_data.population_decoding.(test_name).accuracy = [];
        agg_data.population_decoding.(test_name).accuracy_ci = [];
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

    session_data_path = fullfile(project_root, 'data', 'processed', ...
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
    if isfield(session_data.analysis, 'anova_results')
        for j = 1:length(analysis_plan.anova_plan)
            plan_item = analysis_plan.anova_plan(j);
            analysis_name = plan_item.name;

            if isfield(session_data.analysis.anova_results, analysis_name)
                session_results = session_data.analysis.anova_results.(analysis_name);
                agg_data.anova_results.(analysis_name) = [ ...
                    agg_data.anova_results.(analysis_name); ...
                    session_results ...
                ];
            end
        end
    end

    % --- 4. Aggregate Behavioral Results ---
    if isfield(session_data.analysis, 'behavioral_results')
        for j = 1:length(analysis_plan.behavior_plan)
            plan_item = analysis_plan.behavior_plan(j);
            analysis_name = plan_item.name;

            if isfield(session_data.analysis.behavioral_results, analysis_name)
                session_table = session_data.analysis.behavioral_results.(analysis_name);

                % Add session_id to the table to track provenance
                session_id_col = repmat({session_id}, height(session_table), 1);
                session_table.session_id = session_id_col;

                % Append to the master table
                agg_data.behavioral_results.(analysis_name) = [ ...
                    agg_data.behavioral_results.(analysis_name); ...
                    session_table ...
                ];
            end
        end
    end

    % --- 5. Aggregate Population Decoding Results ---
    if isfield(session_data.analysis, 'population_decoding')
        for j = 1:length(analysis_plan.decoding_plan.testing_plan)
            test = analysis_plan.decoding_plan.testing_plan(j);
            test_name = test.test_name;

            if isfield(session_data.analysis.population_decoding, test_name)
                res = session_data.analysis.population_decoding.(test_name);
                agg_data.population_decoding.(test_name).accuracy = [ ...
                    agg_data.population_decoding.(test_name).accuracy; ...
                    res.accuracy ...
                ];
                agg_data.population_decoding.(test_name).accuracy_ci = [ ...
                    agg_data.population_decoding.(test_name).accuracy_ci; ...
                    res.accuracy_ci ...
                ];
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
