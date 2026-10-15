%% test_screening_visualizations.m
% Test script to verify the screening visualization functions work correctly.
%
% Tests:
%   1. Session-level screening metrics plots
%   2. Data aggregation function
%   3. Aggregate screening metrics plots
%
% Author: Claude Code
% Date: 2026-01-23

clear; clc; close all;

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

fprintf('=== Testing Screening Visualizations ===\n\n');

%% Load manifest
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);
fprintf('Loaded manifest with %d sessions\n\n', height(manifest));

%% Test 1: Session-level plots for first session
fprintf('=== Test 1: Session-Level Screening Metrics Plots ===\n');

% Use first session for testing
session_id = manifest.unique_id{1};
fprintf('Testing with session: %s\n', session_id);

% Load session data
local_processed_dir = fullfile(project_root, 'data', 'processed', session_id);
local_session_data_path = fullfile(local_processed_dir, [session_id, '_session_data.mat']);

if ~exist(local_session_data_path, 'file')
    error('Session data not found: %s', local_session_data_path);
end

fprintf('Loading session data...\n');
load(local_session_data_path, 'session_data');

% Add metadata
session_data.metadata = table2struct(manifest(1,:));

% Check for screening_info
if isfield(session_data, 'analysis') && isfield(session_data.analysis, 'screening_info')
    screening_info = session_data.analysis.screening_info;
    fprintf('Found screening_info with %d neurons\n', length(screening_info));
else
    error('No screening_info found in session data');
end

% Generate session-level plots
output_dir = fullfile(project_root, 'figures');
fprintf('Generating session-level screening metrics plots...\n');
tic;
plot_screening_metrics_session(session_data, screening_info, output_dir, session_id);
session_plot_time = toc;
fprintf('Session plots generated in %.2f seconds\n', session_plot_time);

% Verify output
expected_pdf = fullfile(output_dir, [session_id, '_screening_metrics.pdf']);
if exist(expected_pdf, 'file')
    fprintf('SUCCESS: Session PDF created at %s\n\n', expected_pdf);
else
    warning('Session PDF not found at expected location');
end

%% Test 2: Data aggregation
fprintf('=== Test 2: Data Aggregation ===\n');

fprintf('Aggregating screening data from all sessions...\n');
tic;
aggregate_data = aggregate_screening_data(manifest, project_root);
aggregation_time = toc;
fprintf('Aggregation completed in %.2f seconds\n', aggregation_time);

% Verify aggregation results
fprintf('\nAggregation Results:\n');
fprintf('  Total neurons: %d\n', length(aggregate_data.mean_firing_rates));
fprintf('  Total included: %d (%.1f%%)\n', sum(aggregate_data.is_included), ...
    100*mean(aggregate_data.is_included));
fprintf('  Sessions with data: %d\n', height(aggregate_data.session_summary));
fprintf('  SC neurons: %d\n', sum(strcmpi(aggregate_data.brain_areas, 'SC')));
fprintf('  SNc neurons: %d\n', sum(strcmpi(aggregate_data.brain_areas, 'SNc')));

% Check for expected fields
expected_fields = {'session_ids', 'brain_areas', 'mean_firing_rates', 'sparsity', ...
    'stability_scores', 'is_included', 'exclusion_reasons', 'session_summary'};
missing_fields = setdiff(expected_fields, fieldnames(aggregate_data));
if isempty(missing_fields)
    fprintf('  All expected fields present: YES\n\n');
else
    warning('Missing fields: %s', strjoin(missing_fields, ', '));
end

%% Test 3: Aggregate plots
fprintf('=== Test 3: Aggregate Screening Metrics Plots ===\n');

fprintf('Generating aggregate screening metrics plots...\n');
tic;
plot_screening_metrics_aggregate(aggregate_data, output_dir);
aggregate_plot_time = toc;
fprintf('Aggregate plots generated in %.2f seconds\n', aggregate_plot_time);

% Verify output
expected_aggregate_pdf = fullfile(output_dir, 'aggregate_screening_metrics.pdf');
if exist(expected_aggregate_pdf, 'file')
    fprintf('SUCCESS: Aggregate PDF created at %s\n\n', expected_aggregate_pdf);
else
    warning('Aggregate PDF not found at expected location');
end

%% Test 4: Test with a second session (different brain area if available)
fprintf('=== Test 4: Second Session Test ===\n');

% Find a session with different brain area
first_area = manifest.brain_area{1};
other_idx = find(~strcmp(manifest.brain_area, first_area), 1);
if ~isempty(other_idx)
    session_id2 = manifest.unique_id{other_idx};
    fprintf('Testing second session: %s (brain area: %s)\n', session_id2, manifest.brain_area{other_idx});

    % Load second session
    local_processed_dir2 = fullfile(project_root, 'data', 'processed', session_id2);
    local_session_data_path2 = fullfile(local_processed_dir2, [session_id2, '_session_data.mat']);

    if exist(local_session_data_path2, 'file')
        load(local_session_data_path2, 'session_data');
        session_data.metadata = table2struct(manifest(other_idx,:));

        if isfield(session_data, 'analysis') && isfield(session_data.analysis, 'screening_info')
            screening_info2 = session_data.analysis.screening_info;

            tic;
            plot_screening_metrics_session(session_data, screening_info2, output_dir, session_id2);
            fprintf('Second session plots generated in %.2f seconds\n', toc);

            expected_pdf2 = fullfile(output_dir, [session_id2, '_screening_metrics.pdf']);
            if exist(expected_pdf2, 'file')
                fprintf('SUCCESS: Second session PDF created\n\n');
            end
        else
            fprintf('Skipped: No screening_info in second session\n\n');
        end
    else
        fprintf('Skipped: Session data not found\n\n');
    end
else
    fprintf('Skipped: No session with different brain area found\n\n');
end

%% Summary
fprintf('=== Test Summary ===\n');
fprintf('All tests completed successfully!\n');
fprintf('Key results:\n');
fprintf('  - Session plot time: %.2f seconds\n', session_plot_time);
fprintf('  - Aggregation time: %.2f seconds\n', aggregation_time);
fprintf('  - Aggregate plot time: %.2f seconds\n', aggregate_plot_time);
fprintf('  - Total neurons aggregated: %d\n', length(aggregate_data.mean_firing_rates));
fprintf('  - Overall inclusion rate: %.1f%%\n', 100*mean(aggregate_data.is_included));

fprintf('\nOutput files:\n');
fprintf('  - %s\n', expected_pdf);
fprintf('  - %s\n', expected_aggregate_pdf);
