%% generate_all_screening_plots.m
% Generate screening metrics plots for all sessions.
%
% Author: Claude Code
% Date: 2026-01-23

clear; clc; close all;

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

fprintf('=== Generating Screening Metrics Plots for All Sessions ===\n\n');

%% Load manifest
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);
n_sessions = height(manifest);
fprintf('Found %d sessions in manifest\n\n', n_sessions);

output_dir = fullfile(project_root, 'figures');

%% Generate session-level plots for each session
fprintf('=== Generating Session-Level Plots ===\n');
tic;

for i = 1:n_sessions
    session_id = manifest.unique_id{i};
    fprintf('[%d/%d] Processing %s...\n', i, n_sessions, session_id);

    % Load session data
    local_processed_dir = fullfile(project_root, 'data', 'processed', session_id);
    local_session_data_path = fullfile(local_processed_dir, [session_id, '_session_data.mat']);

    if ~exist(local_session_data_path, 'file')
        fprintf('  Skipping: session data not found\n');
        continue;
    end

    load(local_session_data_path, 'session_data');

    % Add metadata
    session_data.metadata = table2struct(manifest(i,:));

    % Check for screening_info
    if ~isfield(session_data, 'analysis') || ~isfield(session_data.analysis, 'screening_info')
        fprintf('  Skipping: no screening_info found\n');
        continue;
    end

    screening_info = session_data.analysis.screening_info;

    % Generate plots
    plot_screening_metrics_session(session_data, screening_info, output_dir, session_id);
end

session_time = toc;
fprintf('\nSession-level plots completed in %.1f seconds\n\n', session_time);

%% Generate aggregate plots
fprintf('=== Generating Aggregate Plots ===\n');
tic;

aggregate_data = aggregate_screening_data(manifest, project_root);
plot_screening_metrics_aggregate(aggregate_data, output_dir);

aggregate_time = toc;
fprintf('\nAggregate plots completed in %.1f seconds\n\n', aggregate_time);

%% Summary
fprintf('=== Complete ===\n');
fprintf('Total time: %.1f seconds\n', session_time + aggregate_time);
fprintf('Output directory: %s\n', output_dir);
fprintf('\nGenerated files:\n');

% List generated files
pdf_files = dir(fullfile(output_dir, '*screening_metrics*.pdf'));
for i = 1:length(pdf_files)
    fprintf('  - %s\n', pdf_files(i).name);
end
