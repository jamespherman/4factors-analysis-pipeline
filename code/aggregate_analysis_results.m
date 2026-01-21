%% aggregate_analysis_results.m
%
% This script aggregates single-session analysis results into a
% "plot-ready" data structure, `aggregated_analysis_data.mat`. It is
% entirely "plan-driven," meaning it uses the `analysis_plan` from
% `define_task_conditions.m` as the sole source of truth for
% determining the structure of the output.
%
% The script is designed to be robust to session-to-session variability in
% analysis outputs. It standardizes the aggregation process by ensuring
% that for each analysis type, the output struct array has a consistent
% set of fields, even when a given session is missing a particular result.
%
% This is achieved using an "assign-on-first-pass, then-append" logic.
% For the first session of a given brain area, the script directly assigns
% the analysis result to the aggregated data structure. For all subsequent
% sessions, it appends the results. If a session is missing a specific
% analysis, a standardized placeholder struct (or table) is used, which
% guarantees that the final aggregated array is always valid and prevents
% crashes from inconsistent field names.
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

    % 4. Window-Based ROC Results
    if isfield(analysis_plan, 'window_roc_plan')
        agg_data.window_roc = struct();
        template_struct = struct('session_id', '', 'n_neurons', NaN, ...
            'auc', [], 'p', [], 'ci', [], 'n_trials', []);

        epoch_names = fieldnames(analysis_plan.window_roc_plan.epoch_windows);
        factor_names = {analysis_plan.window_roc_plan.factors.name};

        for i_epoch = 1:length(epoch_names)
            epoch_name = epoch_names{i_epoch};
            agg_data.window_roc.(epoch_name) = struct();
            for i_factor = 1:length(factor_names)
                factor_name = factor_names{i_factor};
                agg_data.window_roc.(epoch_name).(factor_name) = template_struct;
            end
        end
    end

    % 5. Neuron Metrics Table (Task 7.2)
    % Initialize an empty table that will be populated with per-neuron data
    % including identifiers, metrics, and window-based ROC values
    agg_data.neuron_metrics_table = table();

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

                % What's the plan name?
                plan_name = plan_item.name;

                % Check if source data exists
                source_exists = ...
                    isfield(analysis_results, 'roc_comparison') && ...
                    isfield(analysis_results.roc_comparison, ...
                    event_name) && ...
                    isfield(analysis_results.roc_comparison ...
                    .(event_name), plan_name);

                % if course data exists put it in 'session_struct'
                if source_exists
                    session_struct = analysis_results.roc_comparison ...
                        .(event_name).(plan_name);
                    session_struct.session_id = session_id;
                    session_struct.n_neurons = n_neurons;
                end

                % If this is the 1st session for this brain area, just
                % assign the data, otherwise append it:
                if length(session_ids.(lower(brain_area))) == 1
                    agg_data.roc_comparison.(event_name).(plan_name) ...
                        = session_struct;
                else
                    agg_data.roc_comparison.(event_name).(plan_name)(...
                        end+1) ...
                        = session_struct;
                end
            end
        end
    end

    % Baseline Comparison
    if isfield(agg_data, 'baseline_comparison')
        for event_cell = analysis_plan.events
            event_name = event_cell{:};
            for plan_item = analysis_plan.baseline_plan

                % What's the plan name?
                plan_name = plan_item.name;

                % Check if source data exists
                source_exists = ...
                    isfield(analysis_results, 'baseline_comparison') && ...
                    isfield(analysis_results.baseline_comparison, ...
                    event_name) && ...
                    isfield(analysis_results.baseline_comparison ...
                    .(event_name), plan_name);

                % if sourse data exists put it in 'session_struct'
                if source_exists
                    session_struct = ...
                        analysis_results.baseline_comparison ...
                        .(event_name).(plan_name);
                    session_struct.session_id = session_id;
                    session_struct.n_neurons = n_neurons;

                end

                % If this is the 1st session for this brain area, just
                % assign the data, otherwise append it:
                if length(session_ids.(lower(brain_area))) == 1
                    agg_data.baseline_comparison.(event_name).(plan_name) ...
                        = session_struct;
                else
                    agg_data.baseline_comparison.(event_name).(plan_name)(end+1) ...
                        = session_struct;
                end
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

                    % If this is the 1st session for this brain area, just
                    % assign the data, otherwise append it:
                    if length(session_ids.(lower(brain_area))) == 1
                        agg_data.behavioral_results.(analysis_name) = ...
                            session_table;
                    else
                        current_table = ...
                            agg_data.behavioral_results.(analysis_name);
                        agg_data.behavioral_results.(analysis_name) = ...
                            [current_table; session_table];
                    end
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

            % Check for source data and populate
            if isfield(analysis_results.population_decoding, test_name)
                session_struct = ...
                    analysis_results.population_decoding.(test_name);
                session_struct.session_id = session_id;
            end

            % If this is the 1st session for this brain area, just
            % assign the data, otherwise append it:
            if length(session_ids.(lower(brain_area))) == 1
                agg_data.population_decoding.(test_name) = ...
                    session_struct;
            else
                agg_data.population_decoding.(test_name)(end+1) = ...
                    session_struct;
            end
        end
    end

    % --- 4. Aggregate Window-Based ROC Results ---
    if isfield(agg_data, 'window_roc') && ...
            isfield(analysis_results, 'window_roc')

        epoch_names = fieldnames(analysis_results.window_roc);
        for i_epoch = 1:length(epoch_names)
            epoch_name = epoch_names{i_epoch};

            if ~isfield(agg_data.window_roc, epoch_name)
                continue;
            end

            factor_names = fieldnames(analysis_results.window_roc.(epoch_name));
            for i_factor = 1:length(factor_names)
                factor_name = factor_names{i_factor};

                if ~isfield(agg_data.window_roc.(epoch_name), factor_name)
                    continue;
                end

                source_data = analysis_results.window_roc.(epoch_name).(factor_name);
                if ~isfield(source_data, 'auc')
                    continue;
                end

                session_struct = source_data;
                session_struct.session_id = session_id;
                session_struct.n_neurons = n_neurons;

                if length(session_ids.(lower(brain_area))) == 1
                    agg_data.window_roc.(epoch_name).(factor_name) = session_struct;
                else
                    agg_data.window_roc.(epoch_name).(factor_name)(end+1) = session_struct;
                end
            end
        end
    end

    % --- 5. Build Neuron Metrics Table (Task 7.2) ---
    % Create per-neuron rows with identifiers, metrics, and ROC values
    giveFeed(sprintf('  Building neuron metrics table for %s...', session_id));

    % Get number of neurons in this session
    n_session_neurons = n_neurons;

    % Extract cluster IDs
    if isfield(data.session_data, 'spikes') && ...
            isfield(data.session_data.spikes, 'cluster_info') && ...
            isfield(data.session_data.spikes.cluster_info, 'cluster_id')
        cluster_ids = data.session_data.spikes.cluster_info.cluster_id;
        if istable(cluster_ids)
            cluster_ids = cluster_ids{:,1};
        end
    else
        cluster_ids = (1:n_session_neurons)';
    end

    % Build identifier columns
    session_id_col = repmat({session_id}, n_session_neurons, 1);
    cluster_id_col = cluster_ids(1:n_session_neurons);
    brain_area_col = repmat({brain_area}, n_session_neurons, 1);

    % Get SNc subregion (empty for SC)
    if isfield(data.session_data, 'metadata') && ...
            isfield(data.session_data.metadata, 'snc_subregion')
        snc_subregion = data.session_data.metadata.snc_subregion;
        if isempty(snc_subregion) || strcmp(brain_area, 'SC')
            snc_subregion_col = repmat({''}, n_session_neurons, 1);
        else
            snc_subregion_col = repmat({snc_subregion}, n_session_neurons, 1);
        end
    else
        snc_subregion_col = repmat({''}, n_session_neurons, 1);
    end

    % Get neuron selection status (neuron_type: selected = true)
    if isfield(analysis_results, 'selected_neurons')
        is_selected = analysis_results.selected_neurons;
        if length(is_selected) < n_session_neurons
            is_selected = [is_selected; false(n_session_neurons - length(is_selected), 1)];
        end
    else
        is_selected = false(n_session_neurons, 1);
    end

    % Get neuron_type from screening_info if available
    % Categories for SC: 'classic', 'interneuron', 'excluded'
    % Categories for SNc: 'modulated', 'excluded'
    % The neuron_class field is populated by apply_neuron_screening
    neuron_type_col = repmat({''}, n_session_neurons, 1);
    if isfield(analysis_results, 'screening_info') && ...
            ~isempty(analysis_results.screening_info)
        screening_info = analysis_results.screening_info;
        for i_neuron = 1:min(length(screening_info), n_session_neurons)
            % Use neuron_class if available (SC-specific classification)
            if isfield(screening_info(i_neuron), 'neuron_class') && ...
                    ~isempty(screening_info(i_neuron).neuron_class)
                neuron_type_col{i_neuron} = screening_info(i_neuron).neuron_class;
            elseif screening_info(i_neuron).included
                neuron_type_col{i_neuron} = 'modulated';
            else
                % Use exclusion_reason if available, otherwise generic 'excluded'
                if isfield(screening_info(i_neuron), 'exclusion_reason') && ...
                        ~isempty(screening_info(i_neuron).exclusion_reason)
                    neuron_type_col{i_neuron} = ['excluded: ' screening_info(i_neuron).exclusion_reason];
                else
                    neuron_type_col{i_neuron} = 'excluded';
                end
            end
        end
    else
        % Fallback: derive from is_selected
        for i_neuron = 1:n_session_neurons
            if is_selected(i_neuron)
                neuron_type_col{i_neuron} = 'modulated';
            else
                neuron_type_col{i_neuron} = 'excluded';
            end
        end
    end

    % Extract additional screening metrics if available
    screening_mean_fr = NaN(n_session_neurons, 1);
    screening_sparsity = NaN(n_session_neurons, 1);
    passes_fr_threshold = false(n_session_neurons, 1);
    passes_sparsity_threshold = false(n_session_neurons, 1);
    is_task_modulated = false(n_session_neurons, 1);

    if isfield(analysis_results, 'screening_info') && ~isempty(analysis_results.screening_info)
        screening_info = analysis_results.screening_info;
        for i_neuron = 1:min(length(screening_info), n_session_neurons)
            if isfield(screening_info(i_neuron), 'mean_firing_rate')
                screening_mean_fr(i_neuron) = screening_info(i_neuron).mean_firing_rate;
            end
            if isfield(screening_info(i_neuron), 'proportion_empty_bins')
                screening_sparsity(i_neuron) = screening_info(i_neuron).proportion_empty_bins;
            end
            if isfield(screening_info(i_neuron), 'passes_fr_threshold')
                passes_fr_threshold(i_neuron) = screening_info(i_neuron).passes_fr_threshold;
            end
            if isfield(screening_info(i_neuron), 'passes_sparsity_threshold')
                passes_sparsity_threshold(i_neuron) = screening_info(i_neuron).passes_sparsity_threshold;
            end
            if isfield(screening_info(i_neuron), 'is_modulated')
                is_task_modulated(i_neuron) = screening_info(i_neuron).is_modulated;
            end
        end
    end

    % Get baseline firing rates
    if isfield(data.session_data, 'metrics') && ...
            isfield(data.session_data.metrics, 'baseline_frs')
        baseline_frs = data.session_data.metrics.baseline_frs(:);
        if length(baseline_frs) < n_session_neurons
            baseline_frs = [baseline_frs; NaN(n_session_neurons - length(baseline_frs), 1)];
        end
    else
        baseline_frs = NaN(n_session_neurons, 1);
    end

    % Get waveform durations (peak-to-trough)
    if isfield(data.session_data, 'metrics') && ...
            isfield(data.session_data.metrics, 'wf_metrics')
        wf_metrics = data.session_data.metrics.wf_metrics;
        if isstruct(wf_metrics)
            wf_durations = [wf_metrics.peak_trough_ms]';
        else
            wf_durations = NaN(n_session_neurons, 1);
        end
        if length(wf_durations) < n_session_neurons
            wf_durations = [wf_durations; NaN(n_session_neurons - length(wf_durations), 1)];
        end
    else
        wf_durations = NaN(n_session_neurons, 1);
    end

    % Build initial table with identifiers and metrics
    session_neuron_table = table(...
        session_id_col, cluster_id_col, brain_area_col, snc_subregion_col, ...
        neuron_type_col, is_selected, baseline_frs, wf_durations, ...
        screening_mean_fr, screening_sparsity, ...
        passes_fr_threshold, passes_sparsity_threshold, is_task_modulated, ...
        'VariableNames', {'session_id', 'cluster_id', 'brain_area', ...
            'snc_subregion', 'neuron_type', 'is_selected', 'baseline_fr', ...
            'waveform_duration', 'screening_mean_fr', 'screening_sparsity', ...
            'passes_fr_threshold', 'passes_sparsity_threshold', 'is_task_modulated'});

    % Add window-based ROC columns for each epoch × factor combination
    if isfield(analysis_results, 'window_roc')
        wroc_epochs = fieldnames(analysis_results.window_roc);
        for i_epoch = 1:length(wroc_epochs)
            epoch_name = wroc_epochs{i_epoch};
            if ~isstruct(analysis_results.window_roc.(epoch_name))
                continue;
            end
            wroc_factors = fieldnames(analysis_results.window_roc.(epoch_name));
            for i_factor = 1:length(wroc_factors)
                factor_name = wroc_factors{i_factor};
                wroc_data = analysis_results.window_roc.(epoch_name).(factor_name);

                % Column names: auc_{epoch}_{factor}, p_{epoch}_{factor}
                auc_col_name = sprintf('auc_%s_%s', epoch_name, factor_name);
                p_col_name = sprintf('p_%s_%s', epoch_name, factor_name);

                % Extract AUC and p-values
                if isfield(wroc_data, 'auc') && ~isempty(wroc_data.auc)
                    auc_vals = wroc_data.auc(:);
                    if length(auc_vals) < n_session_neurons
                        auc_vals = [auc_vals; NaN(n_session_neurons - length(auc_vals), 1)];
                    end
                else
                    auc_vals = NaN(n_session_neurons, 1);
                end

                if isfield(wroc_data, 'p') && ~isempty(wroc_data.p)
                    p_vals = wroc_data.p(:);
                    if length(p_vals) < n_session_neurons
                        p_vals = [p_vals; NaN(n_session_neurons - length(p_vals), 1)];
                    end
                else
                    p_vals = NaN(n_session_neurons, 1);
                end

                % Add columns to table
                session_neuron_table.(auc_col_name) = auc_vals;
                session_neuron_table.(p_col_name) = p_vals;
            end
        end
    end

    % Append to aggregated neuron metrics table
    if isempty(agg_data.neuron_metrics_table)
        agg_data.neuron_metrics_table = session_neuron_table;
    else
        % Ensure tables have same columns before concatenating
        existing_cols = agg_data.neuron_metrics_table.Properties.VariableNames;
        new_cols = session_neuron_table.Properties.VariableNames;

        % Add any missing columns to existing table
        for col_cell = setdiff(new_cols, existing_cols)
            col_name = col_cell{1};
            n_existing = height(agg_data.neuron_metrics_table);
            if isnumeric(session_neuron_table.(col_name))
                agg_data.neuron_metrics_table.(col_name) = NaN(n_existing, 1);
            else
                agg_data.neuron_metrics_table.(col_name) = repmat({''}, n_existing, 1);
            end
        end

        % Add any missing columns to new table
        for col_cell = setdiff(existing_cols, new_cols)
            col_name = col_cell{1};
            if isnumeric(agg_data.neuron_metrics_table.(col_name))
                session_neuron_table.(col_name) = NaN(n_session_neurons, 1);
            else
                session_neuron_table.(col_name) = repmat({''}, n_session_neurons, 1);
            end
        end

        % Reorder columns to match
        session_neuron_table = session_neuron_table(:, agg_data.neuron_metrics_table.Properties.VariableNames);

        % Concatenate
        agg_data.neuron_metrics_table = [agg_data.neuron_metrics_table; session_neuron_table];
    end
    giveFeed(sprintf('  Added %d neurons to metrics table.', n_session_neurons));

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