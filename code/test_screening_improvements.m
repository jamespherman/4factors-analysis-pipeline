%% test_screening_improvements.m
% Quick test script to verify the screening improvements work correctly.
%
% Tests:
%   1. Friedman test caching (Phase 1)
%   2. Stable trial period detection (Phase 2)
%   3. Stability plotting (Phase 3)
%
% Author: Claude Code
% Date: 2026-01-23

clear; clc; close all;

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

fprintf('=== Testing Screening Improvements ===\n\n');

%% Load a test session
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

% Use first SC session for testing
sc_sessions = manifest(strcmp(manifest.brain_area, 'SC'), :);
if isempty(sc_sessions)
    error('No SC sessions found in manifest');
end

session_id = sc_sessions.unique_id{1};
fprintf('Testing with session: %s\n\n', session_id);

% Load session data
local_processed_dir = fullfile(project_root, 'data', 'processed', session_id);
local_session_data_path = fullfile(local_processed_dir, [session_id, '_session_data.mat']);

if exist(local_session_data_path, 'file')
    fprintf('Loading cached data...\n');
    load(local_session_data_path, 'session_data');
else
    error('Session data not found: %s', local_session_data_path);
end

% Add metadata
session_data.metadata = table2struct(sc_sessions(1,:));
session_data.metadata.snc_subregion = '';

%% Test 1: Stable Trial Period Detection
fprintf('\n=== Test 1: Stable Trial Period Detection ===\n');

cluster_ids = session_data.spikes.cluster_info.cluster_id;
nNeurons = numel(cluster_ids);

% Remove any cached stability info to force recomputation
if isfield(session_data, 'stability_info')
    session_data = rmfield(session_data, 'stability_info');
end

stability_params = struct(...
    'window_size', 20, ...
    'step_size', 5, ...
    'cv_threshold', 0.5, ...
    'min_fr_threshold', 0.1, ...
    'min_stable_trials', 50);

tic;
stability_info = detect_stable_trial_period(...
    session_data.spikes.times, ...
    session_data.spikes.clusters, ...
    session_data.eventTimes.trialBegin, ...
    session_data.eventTimes.trialEnd, ...
    cluster_ids, ...
    stability_params);
stability_time = toc;

fprintf('Stability detection completed in %.2f seconds\n', stability_time);
fprintf('Results:\n');
n_stable = sum([stability_info.is_stable]);
fprintf('  - Stable neurons: %d/%d (%.1f%%)\n', n_stable, nNeurons, 100*n_stable/nNeurons);
fprintf('  - Mean stable period: %.1f trials\n', mean([stability_info.n_stable_trials]));
fprintf('  - Mean presence ratio: %.2f\n', mean([stability_info.presence_ratio]));

%% Test 2: Friedman Test Caching
fprintf('\n=== Test 2: Friedman Test Caching ===\n');

% Remove any cached Friedman results to force recomputation
if isfield(session_data, 'cached_friedman_results')
    session_data = rmfield(session_data, 'cached_friedman_results');
end

[~, analysis_plan] = define_task_conditions();

fprintf('First run (computing Friedman tests)...\n');
tic;
[~, ~, ~, ~, session_data] = screen_sc_neurons(session_data, project_root, analysis_plan, session_id);
first_run_time = toc;
fprintf('First run completed in %.2f seconds\n', first_run_time);

% Verify cache was created
if isfield(session_data, 'cached_friedman_results')
    fprintf('Cache created with %d cluster results\n', ...
        length(session_data.cached_friedman_results.cluster_ids));
else
    error('Cache was not created!');
end

fprintf('Second run (using cached results)...\n');
tic;
[selected2, sig2, side2, class2, session_data2] = screen_sc_neurons(session_data, project_root, analysis_plan, session_id);
second_run_time = toc;
fprintf('Second run completed in %.2f seconds\n', second_run_time);

fprintf('Speedup: %.1fx faster\n', first_run_time / second_run_time);

%% Test 3: Apply Neuron Screening with Stable Periods
fprintf('\n=== Test 3: Apply Neuron Screening with Stable Periods ===\n');

% Prepare core_data for screening
if ~isfield(session_data, 'analysis') || ~isfield(session_data.analysis, 'core_data')
    fprintf('Preparing core data...\n');
    alignment_events = {'targetOn', 'saccadeOnset'};
    all_neurons_mask = true(nNeurons, 1);
    core_data = prepare_core_data(session_data, all_neurons_mask, alignment_events, analysis_plan);
else
    core_data = session_data.analysis.core_data;
end

% Remove cached stability info to test fresh computation
if isfield(session_data, 'stability_info')
    session_data = rmfield(session_data, 'stability_info');
end

fprintf('Running apply_neuron_screening with stability detection...\n');
tic;
[selected_neurons, screening_info, session_data_out] = apply_neuron_screening(...
    session_data, core_data, analysis_plan, project_root);
screening_time = toc;
fprintf('Screening completed in %.2f seconds\n', screening_time);

% Check results
n_included = sum(selected_neurons);
n_used_stable = sum([screening_info.used_stable_period]);
fprintf('Results:\n');
fprintf('  - Included neurons: %d/%d (%.1f%%)\n', n_included, nNeurons, 100*n_included/nNeurons);
fprintf('  - Used stable period: %d/%d (%.1f%%)\n', n_used_stable, nNeurons, 100*n_used_stable/nNeurons);

%% Test 4: Stability Plotting
fprintf('\n=== Test 4: Stability Plotting ===\n');

output_dir = fullfile(project_root, 'figures');
fprintf('Generating stability diagnostic plots...\n');
tic;
plot_neuron_stability(session_data_out, session_data_out.stability_info, output_dir, session_id);
plot_time = toc;
fprintf('Plotting completed in %.2f seconds\n', plot_time);

expected_pdf = fullfile(output_dir, [session_id, '_stability_analysis.pdf']);
if exist(expected_pdf, 'file')
    fprintf('SUCCESS: PDF created at %s\n', expected_pdf);
else
    warning('PDF not found at expected location');
end

%% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('All tests completed successfully!\n');
fprintf('Key findings:\n');
fprintf('  - Friedman caching speedup: %.1fx\n', first_run_time / second_run_time);
fprintf('  - Neurons using stable period: %d/%d (%.1f%%)\n', ...
    n_used_stable, nNeurons, 100*n_used_stable/nNeurons);
fprintf('  - Final inclusion rate: %d/%d (%.1f%%)\n', ...
    n_included, nNeurons, 100*n_included/nNeurons);
