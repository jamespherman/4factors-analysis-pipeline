%% aggregate_analysis_results.m
%
% This script aggregates single-session analysis results into a
% "plot-ready" data structure, `aggregated_analysis_data.mat`. It is
% entirely "plan-driven," meaning it uses the `analysis_plan` from
% `define_task_conditions.m` as the sole source of truth for
% determining the structure of the output.
%
% The script is designed to be robust to session-to-session
% variability in analysis outputs. For each analysis type (e.g.,
% ANOVA, ROC), it first creates a "template" struct containing all
% possible data fields initialized to default values (e.g., NaN, []).
% In the main loop, it copies this template for each session,
% populates it with whatever data is available, and then appends the
% standardized struct to the final array. This approach guarantees
% that all structs in the array have the same fields, preventing
% crashes caused by inconsistent struct fields from different
% sessions.
%
% The output file contains four top-level variables:
%   - aggregated_sc_data:  Struct with aggregated data for SC sessions.
%   - aggregated_snc_data: Struct with aggregated data for SNc sessions.
%   - analysis_plan:       A copy of the plan used for aggregation.
%   - session_ids:         Struct with unique session IDs for each area.
%
% The script operates in three main stages:
%   1. Initialization: It defines template structs for each analysis and
%      builds the nested structure of the output variables based on the
%      plan.
%   2. Aggregation Loop: It iterates through each session, copies the
%      relevant template, populates it with the session's analysis
%      results, and appends it to the corresponding struct array.
%   3. Final Save: It saves the populated top-level variables to disk.
%
% Author: Jules
% Date: 2025-09-24 (Refactored to use a template-based approach for
% robust aggregation)
%
function [aggregated_sc_data, aggregated_snc_data] ...
    = aggregate_analysis_results()
% Start timer and define in-line function to give user feedback:
tic;
giveFeed = @(x)disp([num2str(round(toc, 1)) 's - ' x]);

%% Setup Paths and Configuration
giveFeed('Setting up paths and loading configuration...');
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
addpath(project_root);

% Load the session manifest and the canonical analysis plan
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
if ~exist(manifest_path, 'file')
    error('aggregate_analysis_results:manifestNotFound', ...
        'Session manifest not found at: %s', manifest_path);
end
manifest = readtable(manifest_path);
[~, analysis_plan] = define_task_conditions();
giveFeed('Manifest and analysis plan loaded.');

%% Initialize Output Structures
giveFeed('Initializing output data structures based on analysis plan...');
% Initialize the main structs for each brain area
aggregated_sc_data = struct();
aggregated_snc_data = struct();

% Initialize session ID containers
session_ids = struct('sc', {{}}, 'snc', {{}});

% Use a map for easier access during the loop
aggregated_data_map = containers.Map(...
    {'SC', 'SNc'}, {aggregated_sc_data, aggregated_snc_data});

for area_cell = keys(aggregated_data_map)
    area = area_cell{:};
    agg_data = aggregated_data_map(area);

    % 1. Time-Resolved: ANOVA, ROC, Baseline Comparison
    % The target structure is a struct array per comparison.

    % ANOVA Results
    if isfield(analysis_plan, 'anova_plan')
        agg_data.anova_results = struct();
        for plan_item = analysis_plan.anova_plan
            if ~plan_item.run, continue; end
            agg_data.anova_results.(plan_item.name) = struct();
            % Define the template struct for this ANOVA
            template_struct = struct('session_id', '', 'n_neurons', ...
                NaN, 'time_vector', []);
            % Generate all possible interaction terms from factors
            factors = plan_item.factors;
            num_factors = length(factors);
            all_terms = {};
            for k = 1:num_factors
                combs = nchoosek(factors, k);
                for i = 1:size(combs, 1)
                    all_terms{end+1} = strjoin(combs(i,:), '*');
                end
            end
            for term_cell = all_terms
                term_name = matlab.lang.makeValidName(term_cell{:});
                template_struct.(['p_' term_name]) = [];
                template_struct.(['f_' term_name]) = [];
            end
            for event_cell = analysis_plan.events
                event_name = event_cell{:};
                % Initialize the struct array with a template
                agg_data.anova_results.(plan_item.name).(event_name) ...
                    = template_struct;
            end
        end
    end

    % ROC Comparison
    if isfield(analysis_plan, 'roc_plan')
        agg_data.roc_comparison = struct();
        template_struct = struct('sig', [], 'time_vector', [], ...
            'session_id', '', 'n_neurons', NaN);
        for event_cell = analysis_plan.events
            event_name = event_cell{:};
            agg_data.roc_comparison.(event_name) = struct();
            for plan_item = analysis_plan.roc_plan
                agg_data.roc_comparison.(event_name).(plan_item.name) = ...
                    template_struct;
            end
        end
    end

    % Baseline Comparison
    if isfield(analysis_plan, 'baseline_plan')
        agg_data.baseline_comparison = struct();
        template_struct = struct('sig', [], 'time_vector', [], ...
            'session_id', '', 'n_neurons', NaN);
        for event_cell = analysis_plan.events
            event_name = event_cell{:};
            agg_data.baseline_comparison.(event_name) = struct();
            for plan_item = analysis_plan.baseline_plan
                plan_name = plan_item.name;
                agg_data.baseline_comparison.(event_name).(...
                    plan_name) = template_struct;
            end
        end
    end

    % 2. Behavioral Analysis
    % The target structure is a single concatenated table per analysis.
    if isfield(analysis_plan, 'behavior_plan')
        agg_data.behavioral_results = struct();
        for plan_item = analysis_plan.behavior_plan
            agg_data.behavioral_results.(plan_item.name) = table();
        end
    end

    % 3. Population Decoding
    % The target structure is a struct array per test.
    if isfield(analysis_plan, 'decoding_plan') && ...
            isfield(analysis_plan.decoding_plan, 'testing_plan')
        agg_data.population_decoding = struct();
        template_struct = struct('accuracy', NaN, 'accuracy_ci', ...
            [NaN, NaN], 'session_id', '');
        for test_item = analysis_plan.decoding_plan.testing_plan
            test_name = test_item.test_name;
            agg_data.population_decoding.(test_name) = ...
                template_struct;
        end
    end
    aggregated_data_map(area) = agg_data;
end
giveFeed('Data structures initialized.');

%% Main Aggregation Loop
complete_sessions = manifest(strcmp(manifest.analysis_status, ...
    'complete'), :);
nSessions = height(complete_sessions);
giveFeed(sprintf('Starting aggregation for %d completed sessions.', ...
    nSessions));

for i = 1:nSessions
    session_info = complete_sessions(i, :);
    session_id = session_info.unique_id{:};
    brain_area = session_info.brain_area{:};
    giveFeed(sprintf('Processing session %d/%d: %s (%s)', ...
        i, nSessions, session_id, brain_area));

    % Load session data
    session_data_path = fullfile(project_root, 'data', 'processed', ...
        session_id, [session_id '_session_data.mat']);
    if ~exist(session_data_path, 'file')
        warning('aggregate_analysis_results:fileNotFound', ...
            ['Could not find session_data.mat for %s. ' ...
            'Skipping.'], session_id);
        continue;
    end
    data = load(session_data_path, 'session_data');
    if ~isfield(data, 'session_data') || ...
            ~isfield(data.session_data, 'analysis')
        warning('aggregate_analysis_results:noAnalysis', ...
            'No analysis field in data for %s. Skipping.', session_id);
        continue;
    end
    analysis_results = data.session_data.analysis;
    n_neurons = size(analysis_results.selected_neurons, 1);

    % Get the correct aggregated struct for the current brain area
    agg_data = aggregated_data_map(brain_area);

    % Append session ID to the correct list
    if strcmpi(brain_area, 'sc')
        session_ids.sc{end+1} = session_id;
    else
        session_ids.snc{end+1} = session_id;
    end

    % --- 1. Aggregate Time-Resolved Analyses ---
    % ANOVA Results
    if isfield(agg_data, 'anova_results')
        for plan_item = analysis_plan.anova_plan
            if ~plan_item.run, continue; end
            for event_cell = analysis_plan.events

                % What's the current event name?
                event_name = event_cell{:};

                % Check if source data exists
                source_exists = ...
                    isfield(analysis_results, 'anova_results') && ...
                    isfield(analysis_results.anova_results, ...
                    plan_item.name) && ...
                    isfield(analysis_results.anova_results ...
                    .(plan_item.name), event_name);

                % if course data exists put it in 'session_struct'
                if source_exists
                    session_struct = analysis_results.anova_results ...
                        .(plan_item.name).(event_name);

                    try
                        session_struct.session_id = session_id;
                        session_struct.n_neurons = n_neurons;
                    catch me
                        keyboard
                    end
                else
                    % The new part: create a standardized placeholder when data is missing
                    session_struct = create_placeholder_anova_struct(plan_item); % Using a helper function
                    session_struct.session_id = session_id;
                    session_struct.n_neurons = n_neurons; % Still record how many neurons there were
                end

                try
                    % If this is the 1st session for this brain area, just
                    % assign the data, otherwise append it:
                    if length(session_ids.(lower(brain_area))) == 1
                        agg_data.anova_results.(plan_item.name).(...
                            event_name) = session_struct;
                    else
                        agg_data.anova_results.(plan_item.name).(...
                            event_name)(end+1) = session_struct;
                    end

                catch me
                    keyboard
                end
            end
        end
    end

    % ROC Comparison
    if isfield(agg_data, 'roc_comparison')
        for event_cell = analysis_plan.events
            event_name = event_cell{:};
            for plan_item = analysis_plan.roc_plan

                % Create a copy of the template
                session_struct = agg_data.roc_comparison.(event_name) ...
                    .(plan_item.name);
                session_struct.session_id = session_id;
                session_struct.n_neurons = n_neurons;

                % Check for source data and populate
                plan_name = plan_item.name;
                source_exists = ...
                    isfield(analysis_results, 'roc_comparison') && ...
                    isfield(analysis_results.roc_comparison, event_name) && ...
                    isfield(analysis_results.roc_comparison.(event_name), ...
                    plan_name);

                if source_exists
                    source_data = analysis_results.roc_comparison ...
                        .(event_name).(plan_name);
                    session_struct.sig = source_data.sig;
                    session_struct.time_vector = source_data.time_vector;
                end
                current_array = agg_data.roc_comparison.(event_name) ...
                    .(plan_name);
                current_array(end+1) = session_struct;
                agg_data.roc_comparison.(event_name).(plan_name) = ...
                    current_array;
            end
        end
    end

    % Baseline Comparison
    if isfield(agg_data, 'baseline_comparison')
        for event_cell = analysis_plan.events
            event_name = event_cell{:};
            for plan_item = analysis_plan.baseline_plan

                % Create a copy of the template
                session_struct = agg_data.baseline_comparison ...
                    .(event_name).(plan_item.name);
                session_struct.session_id = session_id;
                session_struct.n_neurons = n_neurons;

                % Check for source data and populate
                plan_name = plan_item.name;
                source_exists = ...
                    isfield(analysis_results, 'baseline_comparison') && ...
                    isfield(analysis_results.baseline_comparison, ...
                    event_name) && ...
                    isfield(analysis_results.baseline_comparison ...
                    .(event_name), plan_name);

                if source_exists
                    source_data = analysis_results.baseline_comparison ...
                        .(event_name).(plan_name);
                    session_struct.sig = source_data.sig;
                    session_struct.time_vector = source_data.time_vector;
                end
                current_array = agg_data.baseline_comparison ...
                    .(event_name).(plan_name);
                current_array(end+1) = session_struct;
                agg_data.baseline_comparison.(event_name).(plan_name) = ...
                    current_array;
            end
        end
    end

    % --- 2. Aggregate Behavioral Results ---
    if isfield(agg_data, 'behavioral_results') && ...
            isfield(analysis_results, 'behavioral_results')
        for plan_item = analysis_plan.behavior_plan
            analysis_name = plan_item.name;
            if isfield(analysis_results.behavioral_results, analysis_name)
                session_table = dataset2table(...
                    analysis_results.behavioral_results.(analysis_name));
                if ~isempty(session_table)
                    session_table.session_id = repmat({session_id}, ...
                        height(session_table), 1);

                    current_table = ...
                        agg_data.behavioral_results.(analysis_name);
                    agg_data.behavioral_results.(analysis_name) = ...
                        [current_table; session_table];
                end
            end
        end
    end

    % --- 3. Aggregate Population Decoding Results ---
    if isfield(agg_data, 'population_decoding') && ...
            isfield(analysis_results, 'population_decoding')

        % loop over test items:
        for test_item = analysis_plan.decoding_plan.testing_plan

            % What's the 'test name'?
            test_name = test_item.test_name;

            % Create a copy of the template
            session_struct = ...
                agg_data.population_decoding.(test_name);
            session_struct.session_id = session_id;

            % Check for source data and populate
            if isfield(analysis_results.population_decoding, test_name)
                source_data = ...
                    analysis_results.population_decoding.(test_name);
                session_struct.accuracy = source_data.accuracy;
                session_struct.accuracy_ci = source_data.accuracy_ci;
            end
            current_array = agg_data.population_decoding.(test_name);
            current_array(end+1) = session_struct;
            agg_data.population_decoding.(test_name) = current_array;
        end
    end

    % Update the map with the modified struct
    aggregated_data_map(brain_area) = agg_data;
end
giveFeed('All sessions processed.');

%% Finalize and Save
giveFeed('Saving aggregated data to file...');
% Retrieve the final structs from the map
aggregated_sc_data = aggregated_data_map('SC');
aggregated_snc_data = aggregated_data_map('SNc');

% Define the output file path
save_dir = fullfile(project_root, 'data', 'processed');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save_path = fullfile(save_dir, 'aggregated_analysis_data.mat');

% Save the four top-level variables
save(save_path, ...
    'aggregated_sc_data', ...
    'aggregated_snc_data', ...
    'analysis_plan', ...
    'session_ids', ...
    '-v7.3');

giveFeed(['Aggregation complete. Output saved to: ' save_path]);
end