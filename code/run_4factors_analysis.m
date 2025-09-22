%% run_4factors_analysis.m
%
% Main script to run the entire 4factors task analysis pipeline. It iterates
% through a session manifest, performs neuron screening, prepares core
% data, and runs a series of specified analyses based on a centralized plan.
%
% The script is designed to be idempotent; it checks the status of each
% step in the manifest and skips steps that are already marked 'complete'.
%
% Author: Jules
% Date: 2025-09-13
%

%% Setup
clear; clc; close all;

% --- USER TOGGLES ---
force_rerun = struct(...
    'screening', false, ...
    'diag_pdfs', false, ...
    'dataprep',  false, ...
    'analyses',  false ...
);
% --- END USER TOGGLES ---

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

tic;
giveFeed = @(x)disp([num2str(round(toc, 1)) 's - ' x]);

%% Load Manifest
giveFeed('Loading session manifest...');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
if ~exist(manifest_path, 'file')
    error('run_4factors_analysis:manifestNotFound', ...
        'Session manifest not found at: %s', manifest_path);
end
manifest = readtable(manifest_path);
giveFeed('Manifest loaded.');

%% Load Analysis Plan
giveFeed('Loading analysis plan...');
% The single source of truth for the analysis plan is now
% define_task_conditions. We call it without arguments to get the plan.
[~, analysis_plan] = define_task_conditions();
giveFeed('Analysis plan loaded.');

%% Iterate Through Sessions
for i = 1:height(manifest)
    session_id = manifest.unique_id{i};
    giveFeed(sprintf('--- Starting processing for session: %s ---', ...
        session_id));

    one_drive_path = findOneDrive;
    session_data_path = fullfile(one_drive_path, ...
        'Neuronal Data Analysis', session_id, [session_id ...
        '_session_data.mat']);

    if ~exist(session_data_path, 'file')
        warning('run_4factors_analysis:sessionDataNotFound', ...
            'session_data.mat not found for %s. Skipping.', ...
            session_id);
        continue;
    end

    % tell user data is being loaded, then load data, then tell user data
    % was loaded:
    giveFeed(sprintf('Loading data for %s...', session_id));
    load(session_data_path, 'session_data');
    giveFeed('Data loaded.');

    % make sure we have all our metadata:
    session_data.metadata = table2struct(manifest(i,:));

    % assume that data has not been updated:
    data_updated = false;

    % --- Define Task Conditions for this session ---
    giveFeed('Defining task conditions...');
    conditions = define_task_conditions(session_data);
    giveFeed('Task conditions defined.');

    % --- Dynamically determine a proxy event for checks ---
    all_events = {};
    plan_fields = fieldnames(analysis_plan);
    event_field_names = {'event', 'train_event', 'test_event'};
    for i_field = 1:length(plan_fields)
        sub_plan = analysis_plan.(plan_fields{i_field});
        if isstruct(sub_plan)
            for j = 1:length(sub_plan)
                for k = 1:length(event_field_names)
                    field_name = event_field_names{k};
                    if isfield(sub_plan(j), field_name) && ...
                       ~isempty(sub_plan(j).(field_name))
                        all_events{end+1} = sub_plan(j).(field_name);
                    end
                end
            end
        end
    end
    alignment_events = unique(all_events);
    event_proxy = '';
    if ~isempty(alignment_events)
        event_proxy = alignment_events{1};
    end

    % --- Dry run to calculate total number of steps for this session ---
    n_total_steps = 0;
    unique_id = session_id;
    if ~strcmp(manifest.screening_status{i}, 'complete') || ...
            force_rerun.screening
        n_total_steps = n_total_steps + 1;
    end
    diag_output_dir_dry_run = fullfile(project_root, 'figures', unique_id);
    if (~exist(diag_output_dir_dry_run, 'dir') || ...
            isempty(dir(fullfile(diag_output_dir_dry_run, '*.pdf'))))
        n_total_steps = n_total_steps + 1;
    end
    if ~strcmp(manifest.dataprep_status{i}, 'complete') || ...
            force_rerun.dataprep
        n_total_steps = n_total_steps + 1;
    end

    % --- Dry Run: Count Analysis Steps ---
    giveFeed('Dry run: calculating number of analysis steps...');

    % A. Baseline Comparison Analyses
    if ~isempty(event_proxy)
        for j = 1:length(analysis_plan.baseline_plan)
            comp_name = analysis_plan.baseline_plan(j).name;
            path_to_check = fullfile('analysis', 'baseline_comparison', event_proxy, comp_name);
            S = substruct('.', strsplit(path_to_check, '/'));
            is_missing = false;
            try
                subsref(session_data, S);
            catch
                is_missing = true;
            end
            if is_missing || force_rerun.analyses
                n_total_steps = n_total_steps + 1;
            end
        end
    end

    % B. ROC Comparison Analyses
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);
        path_to_check = fullfile('analysis', 'roc_comparison', comp.event, comp.name);
        S = substruct('.', strsplit(path_to_check, '/'));
        is_missing = false;
        try
            subsref(session_data, S);
        catch
            is_missing = true;
        end
        if is_missing || force_rerun.analyses
            n_total_steps = n_total_steps + 1;
        end
    end

    % C. N-way ANOVA Analysis
    if isfield(analysis_plan, 'anova_plan')
        for j = 1:length(analysis_plan.anova_plan)
            current_plan_item = analysis_plan.anova_plan(j);
            analysis_name = current_plan_item.name;
            path_as_string = ['analysis.anova_results.' analysis_name];
            S = substruct('.', strsplit(path_as_string, '.'));
            is_missing = false;
            try
                subsref(session_data, S);
            catch
                is_missing = true;
            end
            if is_missing || force_rerun.analyses
                n_total_steps = n_total_steps + 1;
            end
        end
    end

    % D. Behavioral Analyses
    % Loop through each behavioral analysis plan.
    if isfield(analysis_plan, 'behavior_plan')
        for j = 1:length(analysis_plan.behavior_plan)
            analysis_name = analysis_plan.behavior_plan(j).name;
            path_to_check = fullfile('analysis', 'behavioral_results', analysis_name);
            S = substruct('.', strsplit(path_to_check, '/'));
            is_missing = false;
            try
                subsref(session_data, S);
            catch
                is_missing = true;
            end
            if is_missing || force_rerun.analyses
                n_total_steps = n_total_steps + 1;
            end
        end
    end

    % E. Population Decoding Analyses
    if isfield(analysis_plan, 'decoding_plan')
        % The execution has two potential steps, which are counted independently.

        % Stage 1: Count training step
        needs_training = force_rerun.analyses || ~isfield(session_data, 'analysis') || ~isfield(session_data.analysis, 'population_decoding');
        if needs_training
            n_total_steps = n_total_steps + 1;
        end

        % Stage 2: Count testing step
        is_any_test_missing = false;
        if isfield(session_data, 'analysis') && isfield(session_data.analysis, 'population_decoding')
            testing_plan = analysis_plan.decoding_plan.testing_plan;
            for j = 1:length(testing_plan)
                test_name = testing_plan(j).test_name;
                if ~isfield(session_data.analysis.population_decoding, test_name)
                    is_any_test_missing = true;
                    break;
                end
            end
        else
            % If the main struct is missing, tests are de facto missing.
            is_any_test_missing = true;
        end

        needs_testing = force_rerun.analyses || is_any_test_missing;
        if needs_testing
            % This step is counted independently of the training step. If
            % training is happening, the whole 'population_decoding' field
            % will be new, so testing will also need to run.
            n_total_steps = n_total_steps + 1;
        end
    end
    step_counter = 0;

    % --- Pipeline Stages ---
    % ... (screening, PDF gen, data prep stages are unchanged) ...
        % --- 1. Neuron Screening ---
    if ~strcmp(manifest.screening_status{i}, 'complete') || ...
        force_rerun.screening
        step_counter = step_counter + 1;
        fprintf(['\n--- Session %s: Starting Step %d of %d: ' ...
            '        Neuron Screening ---\n'], unique_id, ...
            step_counter, n_total_steps);
        giveFeed('Screening status is ''pending''. Running screening...');

        session_data.metadata.unique_id = session_id;
        if strcmp(session_data.metadata.brain_area, 'SNc')
            selected_neurons = screen_da_neurons(session_data, session_id, project_root);
        elseif strcmp(session_data.metadata.brain_area, 'SC')
            [selected_neurons, sig_epoch_comp, scSide] = ...
                screen_sc_neurons(session_data, project_root);
            session_data.analysis.scSide = scSide;
            session_data.analysis.sig_epoch_comparison = sig_epoch_comp;
        else
            warning('run_4factors_analysis:unknownSessionType', ...
                'Unknown brain_area ''%s'' for session %s. Cannot screen neurons.', ...
                session_data.metadata.brain_area, session_id);
            continue;
        end

        session_data.analysis.selected_neurons = selected_neurons;

        data_updated = true; % Mark data as updated
        manifest.screening_status{i} = 'complete';
        giveFeed('Screening complete.');
    else
        giveFeed('Screening already complete. Loading results.');
        selected_neurons = session_data.analysis.selected_neurons;
    end

    % --- Per-Neuron Diagnostic PDF Generation ---
    giveFeed('Checking for per-neuron diagnostic PDFs...');
    diag_output_dir = fullfile(project_root, 'figures', session_id);

    % Check if the directory exists and contains any PDF files
    if (~exist(diag_output_dir, 'dir') || isempty(dir(fullfile( ...
            diag_output_dir, '*.pdf')))) || force_rerun.diag_pdfs
        step_counter = step_counter + 1;
        fprintf(['\n--- Session %s: Starting Step %d of %d: ' ...
            'Diagnostic PDF Generation ---\n'], unique_id, ...
            step_counter, n_total_steps);
        giveFeed('Generating diagnostic PDF...');
        if ~exist(diag_output_dir, 'dir')
            mkdir(diag_output_dir);
        end
        generate_neuron_summary_pdf(session_data, selected_neurons, ...
            session_id, diag_output_dir);
        giveFeed('Diagnostic PDF generation complete.');
    else
        giveFeed('Diagnostic PDF already exists, skipping...');
    end

    % --- 2. Core Data Preparation ---
    if ~strcmp(manifest.dataprep_status{i}, 'complete') || ...
        force_rerun.dataprep
        step_counter = step_counter + 1;
        fprintf(['\n--- Session %s: Starting Step %d of %d: ' ...
            'Core Data Preparation ---\n'], unique_id, ...
            step_counter, n_total_steps);
        giveFeed(['Data prep status is ''pending''. ' ...
            '            Running prepare_core_data...']);

        % --- Integration of Analysis Plan for Data Prep ---
        % Call define_task_conditions without arguments to get the plan.
        [~, analysis_plan_for_prep] = define_task_conditions();

        % Dynamically generate the list of required alignment events from the plan.
        all_events_prep = {};
        plan_fields_prep = fieldnames(analysis_plan_for_prep);
        event_field_names_prep = {'event', 'train_event', 'test_event'};
        for i_field = 1:length(plan_fields_prep)
            sub_plan = analysis_plan_for_prep.(plan_fields_prep{i_field});
            if isstruct(sub_plan)
                for j = 1:length(sub_plan)
                    for k = 1:length(event_field_names_prep)
                        field_name = event_field_names_prep{k};
                        if isfield(sub_plan(j), field_name) && ...
                           ~isempty(sub_plan(j).(field_name))
                            all_events_prep{end+1} = sub_plan(j).(field_name);
                        end
                    end
                end
            end
        end
        alignment_events_prep = unique(all_events_prep);

        % Now call prepare_core_data with the dynamic events.
        core_data = prepare_core_data(session_data, selected_neurons, alignment_events_prep);
        session_data.analysis.core_data = core_data;

        data_updated = true; % Mark data as updated
        manifest.dataprep_status{i} = 'complete';
        giveFeed('Data prep complete.');
    else
        giveFeed('Data prep already complete. Loading core_data.');
        core_data = session_data.analysis.core_data;
    end

    % --- On-Demand Analysis Execution ---
    giveFeed('Checking for missing analyses...');

    % A. Baseline Comparison Analyses
    if ~isempty(event_proxy)
        for j = 1:length(analysis_plan.baseline_plan)
            comp_name = analysis_plan.baseline_plan(j).name;

            % Standardized Idempotency Check
            is_missing = false;
            path_to_check = fullfile('analysis', 'baseline_comparison', event_proxy, comp_name);
            S = substruct('.', strsplit(path_to_check, '/'));
            try
                subsref(session_data, S);
            catch
                is_missing = true;
            end

            if is_missing || force_rerun.analyses
                step_counter = step_counter + 1;
                fprintf(['\n--- Session %s: Step %d/%d: Baseline ' ...
                    'Comparison for %s ---\n'], unique_id, step_counter, ...
                    n_total_steps, comp_name);
                giveFeed(sprintf('--> Running Baseline Comparison: %s', comp_name));

                result_by_event = analyze_baseline_comparison(...
                    core_data, conditions, 'condition', comp_name);

                event_names = fieldnames(result_by_event);
                for k = 1:length(event_names)
                    event_name = event_names{k};
                    session_data.analysis.baseline_comparison.(event_name).(comp_name) = ...
                        result_by_event.(event_name).(comp_name);
                end
                data_updated = true;
            end
        end
    end

    % B. ROC Comparison Analyses
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);

        % Standardized Idempotency Check
        is_missing = false;
        path_to_check = fullfile('analysis', 'roc_comparison', comp.event, comp.name);
        S = substruct('.', strsplit(path_to_check, '/'));
        try
            subsref(session_data, S);
        catch
            is_missing = true;
        end

        if is_missing || force_rerun.analyses
            step_counter = step_counter + 1;
            fprintf(['\n--- Session %s: Step %d/%d: ROC Comparison ' ...
                'for %s ---\n'], unique_id, step_counter, ...
                n_total_steps, comp.name);
            giveFeed(sprintf('--> Running ROC Comparison: %s', comp.name));

            result = analyze_roc_comparison(core_data, conditions, 'comparison', comp);
            session_data.analysis.roc_comparison.(comp.event).(comp.name) = result;
            data_updated = true;
        end
    end

    % C. N-way ANOVA Analysis
    if isfield(analysis_plan, 'anova_plan')
        for j = 1:length(analysis_plan.anova_plan)
            current_plan_item = analysis_plan.anova_plan(j);
            analysis_name = current_plan_item.name;

            % Standardized Idempotency Check
            is_missing = false;
            path_as_string = ['analysis.anova_results.' analysis_name];
            S = substruct('.', strsplit(path_as_string, '.'));
            try
                subsref(session_data, S);
            catch
                is_missing = true;
            end

            if is_missing || force_rerun.analyses
                step_counter = step_counter + 1;
                fprintf('\n--- Session %s: Step %d/%d: N-way ANOVA for %s ---\n', ...
                    unique_id, step_counter, n_total_steps, analysis_name);
                giveFeed(sprintf('--> Running N-way ANOVA: %s', analysis_name));
                % Pass the specific plan item to the analysis function
                session_data = analyze_anova(session_data, core_data, ...
                    conditions, current_plan_item);
                data_updated = true;
            end
        end
    end

    % D. Behavioral Analyses
    % Loop through each behavioral analysis plan.
    if isfield(analysis_plan, 'behavior_plan')
        for j = 1:length(analysis_plan.behavior_plan)
            plan_item = analysis_plan.behavior_plan(j);
            analysis_name = plan_item.name;

            % Standardized Idempotency Check
            is_missing = false;
            path_to_check = fullfile('analysis', 'behavioral_results', analysis_name);
            S = substruct('.', strsplit(path_to_check, '/'));
            try
                subsref(session_data, S);
            catch
                is_missing = true;
            end

            if is_missing || force_rerun.analyses
                step_counter = step_counter + 1;
                fprintf(['\n--- Session %s: Step %d/%d: Behavioral ' ...
                    'Analysis for %s ---\n'], unique_id, step_counter, ...
                    n_total_steps, analysis_name);
                giveFeed(sprintf('--> Running Behavioral Analysis: %s', analysis_name));

                result = analyze_behavior(session_data, conditions, plan_item);
                session_data.analysis.behavioral_results.(analysis_name) = result;
                data_updated = true;
            end
        end
    end

    % E. Population Decoding Analyses
    if isfield(analysis_plan, 'decoding_plan')
        % Idempotency check for the entire decoding suite. The suite is run
        % if any of its final test results are missing, ensuring that
        % transient trained models are always available for testing.
        run_decoding_suite = false;
        if force_rerun.analyses
            run_decoding_suite = true;
        else
            if ~isfield(session_data, 'analysis') || ~isfield(session_data.analysis, 'population_decoding')
                run_decoding_suite = true;
            else
                testing_plan = analysis_plan.decoding_plan.testing_plan;
                for j = 1:length(testing_plan)
                    test_name = testing_plan(j).test_name;
                    if ~isfield(session_data.analysis.population_decoding, test_name)
                        run_decoding_suite = true;
                        break;
                    end
                end
            end
        end

        if run_decoding_suite
            % --- Stage 1: Train Models ---
            step_counter = step_counter + 1;
            fprintf(['\n--- Session %s: Step %d/%d: Population ' ...
                'Decoding: Model Training ---\n'], unique_id, ...
                step_counter, n_total_steps);
            giveFeed('--> Training all decoding models...');

            trained_models = {};
            training_plan = analysis_plan.decoding_plan.training_plan;
            for j = 1:length(training_plan)
                modelInfo = train_decoder(session_data, conditions, ...
                    core_data, training_plan(j));
                trained_models{end+1} = modelInfo;
            end
            giveFeed('--> Model training complete.');

            % --- Stage 2: Test Models ---
            step_counter = step_counter + 1;
            fprintf(['\n--- Session %s: Step %d/%d: Population ' ...
                'Decoding: Model Testing ---\n'], unique_id, ...
                step_counter, n_total_steps);
            giveFeed('--> Running decoding model tests...');

            testing_plan = analysis_plan.decoding_plan.testing_plan;
            for j = 1:length(testing_plan)
                testing_item = testing_plan(j);
                results = test_decoder(trained_models, conditions, core_data, testing_item);
                session_data.analysis.population_decoding.(testing_item.test_name) = results;
                data_updated = true;
            end
            giveFeed('--> Model testing complete.');
        else
            giveFeed('--> Population decoding suite already complete.');
        end
    end

    % --- Save Updated Data ---
    if data_updated
        giveFeed('Data was updated, saving back to session_data.mat...');
        save(session_data_path, 'session_data', '-v7.3');
        giveFeed('Save complete.');
    else
        giveFeed(['No new analyses or processing were required for ' ...
            '                this session.']);
    end

    % --- Verify Analysis Completion & Update Manifest ---
    giveFeed('Verifying analysis completion status...');
    is_analysis_complete = true;

    % A. Check baseline comparisons
    if is_analysis_complete && ~isempty(event_proxy)
        for j = 1:length(analysis_plan.baseline_plan)
            comp_name = analysis_plan.baseline_plan(j).name;
            path_to_check = fullfile('analysis', 'baseline_comparison', event_proxy, comp_name);
            S = substruct('.', strsplit(path_to_check, '/'));
            try
                subsref(session_data, S);
            catch
                is_analysis_complete = false; break;
            end
        end
    end

    % B. Check ROC comparisons
    if is_analysis_complete
        for j = 1:length(analysis_plan.roc_plan)
            comp = analysis_plan.roc_plan(j);
            path_to_check = fullfile('analysis', 'roc_comparison', comp.event, comp.name);
            S = substruct('.', strsplit(path_to_check, '/'));
            try
                subsref(session_data, S);
            catch
                is_analysis_complete = false; break;
            end
        end
    end

    % C. Check N-way ANOVA
    if is_analysis_complete && isfield(analysis_plan, 'anova_plan')
        for j = 1:length(analysis_plan.anova_plan)
            current_plan_item = analysis_plan.anova_plan(j);
            analysis_name = current_plan_item.name;
            path_as_string = ['analysis.anova_results.' analysis_name];
            S = substruct('.', strsplit(path_as_string, '.'));
            try
                subsref(session_data, S);
            catch
                is_analysis_complete = false;
                break;
            end
        end
    end

    % D. Check behavioral analyses
    % Loop through each behavioral analysis plan to verify completion.
    if is_analysis_complete && isfield(analysis_plan, 'behavior_plan')
        for j = 1:length(analysis_plan.behavior_plan)
            analysis_name = analysis_plan.behavior_plan(j).name;
            path_to_check = fullfile('analysis', 'behavioral_results', analysis_name);
            S = substruct('.', strsplit(path_to_check, '/'));
            try
                subsref(session_data, S);
            catch
                is_analysis_complete = false; break;
            end
        end
    end

    % E. Check population decoding results
    if is_analysis_complete && isfield(analysis_plan, 'decoding_plan')
        testing_plan = analysis_plan.decoding_plan.testing_plan;
        for j = 1:length(testing_plan)
            test_name = testing_plan(j).test_name;
            path_to_check = fullfile('analysis', 'population_decoding', test_name);
            S = substruct('.', strsplit(path_to_check, '/'));
            try
                subsref(session_data, S);
            catch
                is_analysis_complete = false; break;
            end
        end
    end

    if is_analysis_complete
        manifest.analysis_status{i} = 'complete';
        giveFeed('All analyses for this session are complete.');
    else
        giveFeed('Some analyses for this session are still pending.');
    end

    giveFeed(sprintf('--- Finished processing for session: %s ---\n', ...
        session_id));
end

%% Finalize
giveFeed('All sessions processed. Saving updated manifest...');
writetable(manifest, manifest_path);
giveFeed('Manifest saved. Total time elapsed:');
toc
