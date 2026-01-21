%% test_screening_implementation.m
%
% Test script to verify the refined neuron screening implementation.
% Tests on a single session to verify:
% 1. compute_task_modulation works correctly
% 2. apply_neuron_screening applies both criteria
% 3. screening_info is properly populated
% 4. generate_neuron_summary_pdf displays correct headers
%
% Author: Claude Code
% Date: 2026-01-18

clear; clc; close all;

%% Setup paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

%% Load a test session
% Choose one SC session and one SNc session for testing
test_sessions = {'Feynman_08_12_2025_SC', 'Feynman_08_12_2025_SNc'};

for i_session = 1:length(test_sessions)
    session_id = test_sessions{i_session};
    fprintf('\n========================================\n');
    fprintf('Testing session: %s\n', session_id);
    fprintf('========================================\n');

    % Load session data
    local_processed_dir = fullfile(project_root, 'data', 'processed', session_id);
    local_session_data_path = fullfile(local_processed_dir, [session_id, '_session_data.mat']);

    if ~exist(local_session_data_path, 'file')
        % Try OneDrive location
        one_drive_path = findOneDrive;
        source_path = fullfile(one_drive_path, 'Neuronal Data Analysis', ...
            session_id, [session_id '_session_data.mat']);
        if exist(source_path, 'file')
            fprintf('Loading from OneDrive: %s\n', source_path);
            load(source_path, 'session_data');
        else
            warning('Session data not found for %s. Skipping.', session_id);
            continue;
        end
    else
        fprintf('Loading from local cache: %s\n', local_session_data_path);
        load(local_session_data_path, 'session_data');
    end

    % Load manifest to get metadata
    manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
    manifest = readtable(manifest_path);
    row_idx = find(strcmp(manifest.unique_id, session_id));
    if isempty(row_idx)
        warning('Session %s not found in manifest. Skipping.', session_id);
        continue;
    end
    session_data.metadata = table2struct(manifest(row_idx, :));
    session_data.metadata.unique_id = session_id;

    % Get analysis plan
    [~, analysis_plan] = define_task_conditions();

    % Verify metrics exist
    if ~isfield(session_data, 'metrics') || ~isfield(session_data.metrics, 'baseline_frs')
        warning('Session %s missing baseline_frs. Computing on-the-fly...', session_id);
        codes = initCodes();
        config = analysis_plan.baseline_fr_config;
        task_code = codes.(['uniqueTaskCode_' config.task_name]);
        trial_mask = session_data.trialInfo.taskCode == task_code;
        baseline_fr = calculate_baseline_fr(session_data, config.event, config.window, trial_mask);
        session_data.metrics.baseline_frs = baseline_fr;
    end

    %% Step 1: Prepare core_data for ALL neurons
    fprintf('\n--- Step 1: Preparing core_data for ALL neurons ---\n');
    nClusters = height(session_data.spikes.cluster_info);
    all_neurons_mask = true(nClusters, 1);
    alignment_events = analysis_plan.events;

    fprintf('Number of clusters: %d\n', nClusters);
    fprintf('Alignment events: %s\n', strjoin(alignment_events, ', '));

    tic;
    core_data = prepare_core_data(session_data, all_neurons_mask, alignment_events, analysis_plan);
    fprintf('Core data prepared in %.2f seconds.\n', toc);

    % Verify core_data structure
    assert(isfield(core_data, 'spikes'), 'core_data missing spikes field');
    assert(isfield(core_data.spikes, 'targetOn'), 'core_data.spikes missing targetOn');
    fprintf('core_data.spikes.targetOn.rates size: %s\n', mat2str(size(core_data.spikes.targetOn.rates)));

    %% Step 2: Test compute_task_modulation
    fprintf('\n--- Step 2: Testing compute_task_modulation ---\n');
    tic;
    modulation_info = compute_task_modulation(session_data, core_data, analysis_plan);
    fprintf('Task modulation computed in %.2f seconds.\n', toc);

    % Display results
    n_modulated = sum([modulation_info.is_modulated]);
    fprintf('Modulated neurons: %d/%d (%.1f%%)\n', n_modulated, nClusters, 100*n_modulated/nClusters);

    % Show some examples of modulated neurons
    mod_indices = find([modulation_info.is_modulated]);
    if ~isempty(mod_indices)
        fprintf('\nExamples of modulated neurons:\n');
        for j = 1:min(3, length(mod_indices))
            idx = mod_indices(j);
            fprintf('  Neuron %d: %s\n', modulation_info(idx).cluster_id, ...
                strjoin(modulation_info(idx).modulated_at, ', '));
        end
    end

    %% Step 3: Test apply_neuron_screening
    fprintf('\n--- Step 3: Testing apply_neuron_screening ---\n');
    tic;
    [selected_neurons, screening_info] = apply_neuron_screening(session_data, core_data, analysis_plan);
    fprintf('Screening completed in %.2f seconds.\n', toc);

    % Display results
    n_included = sum(selected_neurons);
    fprintf('\nScreening Results:\n');
    fprintf('  Included: %d/%d (%.1f%%)\n', n_included, nClusters, 100*n_included/nClusters);

    % Count by exclusion reason
    n_low_fr = sum(strcmp({screening_info.exclusion_reason}, 'low_fr'));
    n_no_mod = sum(strcmp({screening_info.exclusion_reason}, 'no_mod'));
    n_both = sum(strcmp({screening_info.exclusion_reason}, 'low_fr_and_no_mod'));
    fprintf('  Excluded - low FR only: %d\n', n_low_fr);
    fprintf('  Excluded - no modulation only: %d\n', n_no_mod);
    fprintf('  Excluded - both: %d\n', n_both);

    % Show firing rates for context
    brain_area = session_data.metadata.brain_area;
    if strcmp(brain_area, 'SC')
        fr_threshold = analysis_plan.neuron_inclusion.min_firing_rate_sc;
    else
        fr_threshold = analysis_plan.neuron_inclusion.min_firing_rate_snc;
    end
    fprintf('\nFiring rate distribution (threshold = %.1f sp/s for %s):\n', fr_threshold, brain_area);
    baseline_frs = session_data.metrics.baseline_frs;
    fprintf('  Min: %.2f, Max: %.2f, Median: %.2f sp/s\n', ...
        min(baseline_frs), max(baseline_frs), median(baseline_frs));
    fprintf('  Pass FR threshold: %d/%d\n', sum(baseline_frs >= fr_threshold), nClusters);

    %% Step 4: Test PDF generation (single neuron)
    fprintf('\n--- Step 4: Testing PDF generation ---\n');
    test_output_dir = fullfile(project_root, 'figures', 'test_screening');
    if ~exist(test_output_dir, 'dir')
        mkdir(test_output_dir);
    end

    % Generate PDF for first neuron only (to save time)
    fprintf('Generating test PDF for cluster index 1...\n');
    tic;
    generate_neuron_summary_pdf(session_data, selected_neurons, session_id, ...
        test_output_dir, 'ClusterIndex', 1, 'ScreeningInfo', screening_info);
    fprintf('PDF generated in %.2f seconds.\n', toc);

    % Check if PDF was created
    cluster_id = session_data.spikes.cluster_info.cluster_id(1);
    expected_pdf = fullfile(test_output_dir, sprintf('%s_%03d_summary.pdf', session_id, cluster_id));
    if exist(expected_pdf, 'file')
        fprintf('SUCCESS: PDF created at %s\n', expected_pdf);
    else
        warning('PDF not found at expected location: %s', expected_pdf);
    end

    fprintf('\n');
end

%% Summary
fprintf('\n========================================\n');
fprintf('TEST COMPLETE\n');
fprintf('========================================\n');
fprintf('Please manually inspect the generated PDFs in:\n');
fprintf('  %s\n', test_output_dir);
fprintf('Verify that:\n');
fprintf('  1. Header shows INCLUDED/EXCLUDED with correct status\n');
fprintf('  2. Modulation annotation appears below status\n');
fprintf('  3. Colors are correct (green=included, red=excluded)\n');
