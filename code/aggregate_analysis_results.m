%% aggregate_analysis_results.m
%
% Gathers all single-session analysis results and pools them by
% brain area. This script is fully plan-driven, using the canonical
% analysis_plan from define_task_conditions.m as the single source of
% truth for what data to expect and how to aggregate it.
%
% OUTPUT:
%   aggregated_sc_data  - A struct with aggregated data for SC. Its
%                         structure is a mirror of the analysis_plan.
%   aggregated_snc_data - A struct with aggregated data for SNc. Its
%                         structure is a mirror of the analysis_plan.
%
% Author: Jules
% Date: 2025-09-23 (Refactored to be fully plan-driven)
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
giveFeed('Manifest and plan loaded.');

%% Initialize Aggregated Data Structures
giveFeed('Initializing data structures based on analysis plan...');

brain_areas = {'SC', 'SNc'};
aggregated_data_by_area = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i_area = 1:length(brain_areas)
    area = brain_areas{i_area};
    agg_data = struct();

    % This field will be populated with neuron-level session IDs
    agg_data.session_id = {};

    % 1. Initialize ROC Comparison Results
    % These are expected to be nested by event, then by comparison name.
    if isfield(analysis_plan, 'roc_plan')
        agg_data.roc_comparison = struct();
        for j = 1:length(analysis_plan.roc_plan)
            plan_item = analysis_plan.roc_plan(j);
            % NOTE: roc_plan still has a hardcoded .event field.
            % We will respect that until the plan itself is changed.
            agg_data.roc_comparison.(plan_item.event).(plan_item.name).sig = [];
        end
    end

    % 2. Initialize Baseline Comparison Results
    % These are run for every event in the master event list.
    if isfield(analysis_plan, 'baseline_plan')
        agg_data.baseline_comparison = struct();
        for j = 1:length(analysis_plan.baseline_plan)
            plan_item = analysis_plan.baseline_plan(j);
            for k = 1:length(analysis_plan.events)
                event_name = analysis_plan.events{k};
                agg_data.baseline_comparison.(event_name).(plan_item.name).sig = [];
            end
        end
    end

    % 3. Initialize ANOVA Results
    % These are run for every event in the master event list.
    if isfield(analysis_plan, 'anova_plan')
        agg_data.anova_results = struct();
        for j = 1:length(analysis_plan.anova_plan)
            plan_item = analysis_plan.anova_plan(j);
            if plan_item.run
                agg_data.anova_results.(plan_item.name) = struct();
                for k = 1:length(analysis_plan.events)
                    event_name = analysis_plan.events{k};
                    % Initialize as an empty struct to hold table and stats
                    agg_data.anova_results.(plan_item.name).(event_name) = struct();
                end
            end
        end
    end

    % 4. Initialize Behavioral Results
    if isfield(analysis_plan, 'behavior_plan')
        agg_data.behavioral_results = struct();
        for j = 1:length(analysis_plan.behavior_plan)
            plan_item = analysis_plan.behavior_plan(j);
            agg_data.behavioral_results.(plan_item.name) = table();
        end
    end

    % 5. Initialize Population Decoding Results
    if isfield(analysis_plan, 'decoding_plan') && isfield(analysis_plan.decoding_plan, 'testing_plan')
        agg_data.population_decoding = struct();
        for j = 1:length(analysis_plan.decoding_plan.testing_plan)
            test_name = analysis_plan.decoding_plan.testing_plan(j).test_name;
            agg_data.population_decoding.(test_name).accuracy = [];
            agg_data.population_decoding.(test_name).accuracy_ci = [];
        end
    end

    aggregated_data_by_area(area) = agg_data;
end
giveFeed('Data structures initialized.');


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

    analysis_results = session_data.analysis;
    n_neurons = size(analysis_results.selected_neurons, 1);

    % Get the correct aggregated struct for the current brain area
    agg_data = aggregated_data_by_area(brain_area);

    % Append session ID, one entry per neuron
    agg_data.session_id = [agg_data.session_id; repmat({session_id}, n_neurons, 1)];

    try
    % --- 1. Aggregate ROC Comparison Results ---
    if isfield(analysis_plan, 'roc_plan')
        for j = 1:length(analysis_plan.roc_plan)
            plan_item = analysis_plan.roc_plan(j);
            data_to_append = nan(n_neurons, 1); % Default to NaN
            if isfield(analysis_results, 'roc_comparison') && ...
                    isfield(analysis_results.roc_comparison, plan_item.event) && ...
                    isfield(analysis_results.roc_comparison.(plan_item.event), plan_item.name)

                res = analysis_results.roc_comparison.(plan_item.event).(plan_item.name);
                data_to_append = res.sig;

                % Store time vector if it's the first time we see it
                if ~isfield(agg_data.roc_comparison.(plan_item.event).(plan_item.name), 'time_vector')
                    agg_data.roc_comparison.(plan_item.event).(plan_item.name).time_vector = res.time_vector;
                end
            end
            agg_data.roc_comparison.(plan_item.event).(plan_item.name).sig = ...
                [agg_data.roc_comparison.(plan_item.event).(plan_item.name).sig; data_to_append];
        end
    end
    catch me
        keyboard
    end

    % --- 2. Aggregate Baseline Comparison Results ---
    if isfield(analysis_plan, 'baseline_plan')
        for j = 1:length(analysis_plan.baseline_plan)
            plan_item = analysis_plan.baseline_plan(j);
            for k = 1:length(analysis_plan.events)
                event_name = analysis_plan.events{k};
                data_to_append = nan(n_neurons, 1); % Default to NaN
                if isfield(analysis_results, 'baseline_comparison') && ...
                        isfield(analysis_results.baseline_comparison, event_name) && ...
                        isfield(analysis_results.baseline_comparison.(event_name), plan_item.name)

                    res = analysis_results.baseline_comparison.(event_name).(plan_item.name);
                    data_to_append = res.sig;

                    if ~isfield(agg_data.baseline_comparison.(event_name).(plan_item.name), 'time_vector')
                        agg_data.baseline_comparison.(event_name).(plan_item.name).time_vector = res.time_vector;
                    end
                end
                agg_data.baseline_comparison.(event_name).(plan_item.name).sig = ...
                    [agg_data.baseline_comparison.(event_name).(plan_item.name).sig; data_to_append];
            end
        end
    end

    % --- 3. Aggregate ANOVA Results ---
    if isfield(analysis_plan, 'anova_plan')
        for j = 1:length(analysis_plan.anova_plan)
            plan_item = analysis_plan.anova_plan(j);
            if ~plan_item.run, continue; end

            for k = 1:length(analysis_plan.events)
                event_name = analysis_plan.events{k};
                % For ANOVA, we append the entire result struct. If it's
                % missing, we must append an empty struct to maintain the
                % integrity of the struct array for concatenation.
                data_to_append = struct(); % Default to empty
                if isfield(analysis_results, 'anova_results') && ...
                        isfield(analysis_results.anova_results, plan_item.name) && ...
                        isfield(analysis_results.anova_results.(plan_item.name), event_name)

                    data_to_append = analysis_results.anova_results.(plan_item.name).(event_name);
                end

                % If the aggregated field is currently empty and we are adding
                % the first real result, we can directly assign it. Otherwise,
                % we concatenate.
                if isempty(fieldnames(agg_data.anova_results.(plan_item.name).(event_name)))
                    agg_data.anova_results.(plan_item.name).(event_name) = data_to_append;
                else
                    agg_data.anova_results.(plan_item.name).(event_name) = ...
                        [agg_data.anova_results.(plan_item.name).(event_name); data_to_append];
                end
            end
        end
    end
    try
        % --- 4. Aggregate Behavioral Results ---
        if isfield(analysis_plan, 'behavior_plan') && isfield(analysis_results, 'behavioral_results')
            for j = 1:length(analysis_plan.behavior_plan)
                plan_item = analysis_plan.behavior_plan(j);
                analysis_name = plan_item.name;

                if isfield(analysis_results.behavioral_results, analysis_name)
                    session_table = analysis_results.behavioral_results.(analysis_name);
                    session_id_col = repmat({session_id}, height(session_table), 1);
                    session_table.session_id = session_id_col;

                    % Check if the aggregated table is still empty
                    if isempty(agg_data.behavioral_results.(analysis_name))
                        % If it's the first table for this analysis, just assign it
                        agg_data.behavioral_results.(analysis_name) = session_table;
                    else
                        % For all subsequent tables, perform the vertical concatenation
                        agg_data.behavioral_results.(analysis_name) = [ ...
                            agg_data.behavioral_results.(analysis_name); ...
                            session_table ...
                            ];
                    end
                end
            end
        end
    catch me
        keyboard
    end

    % --- 5. Aggregate Population Decoding Results ---
    if isfield(analysis_plan, 'decoding_plan') && isfield(analysis_results, 'population_decoding')
        for j = 1:length(analysis_plan.decoding_plan.testing_plan)
            test_name = analysis_plan.decoding_plan.testing_plan(j).test_name;

            if isfield(analysis_results.population_decoding, test_name)
                res = analysis_results.population_decoding.(test_name);
                agg_data.population_decoding.(test_name).accuracy = [ ...
                    agg_data.population_decoding.(test_name).accuracy; ...
                    res.accuracy ...
                    ];
                agg_data.population_decoding.(test_name).accuracy_ci = [ ...
                    agg_data.population_decoding.(test_name).accuracy_ci; ...
                    res.accuracy_ci ...
                    ];
            else
                % If a session is missing a decoding result, append NaN
                agg_data.population_decoding.(test_name).accuracy = [ ...
                    agg_data.population_decoding.(test_name).accuracy; nan];
                agg_data.population_decoding.(test_name).accuracy_ci = [ ...
                    agg_data.population_decoding.(test_name).accuracy_ci; nan(1,2)];
            end
        end
    end

    % Update the map with the modified struct
    aggregated_data_by_area(brain_area) = agg_data;
end

giveFeed('All sessions processed.');

%% Finalize and Save
giveFeed('Saving aggregated data to file...');
aggregated_sc_data = aggregated_data_by_area('SC');
aggregated_snc_data = aggregated_data_by_area('SNc');

saveFileName = fullfile(project_root, 'data', 'processed', ...
    'aggregated_analysis_data.mat');
save(saveFileName, 'aggregated_sc_data', 'aggregated_snc_data', '-v7.3');
giveFeed('Aggregation complete.');

end
