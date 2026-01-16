%% test_single_neuron_pdf.m
% Quick test script to generate PDFs for a single session

clear; clc; close all;

% Setup paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

%% Configuration
TEST_SESSION = 'Newton_08_13_2025_SC';  % Change this to test different sessions
MAX_NEURONS = 5;  % Only process first N neurons for quick test

%% Load session data
fprintf('Loading session: %s\n', TEST_SESSION);

session_data_path = fullfile(project_root, 'data', 'processed', ...
    TEST_SESSION, [TEST_SESSION '_session_data.mat']);

if ~exist(session_data_path, 'file')
    error('Session data not found: %s', session_data_path);
end

load(session_data_path, 'session_data');
fprintf('Loaded session data.\n');

%% Load manifest to get metadata
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);
session_row = strcmp(manifest.unique_id, TEST_SESSION);
session_data.metadata = table2struct(manifest(session_row, :));

% Add SNc subregion if applicable
if strcmp(session_data.metadata.brain_area, 'SNc')
    session_data.metadata.snc_subregion = parse_snc_subregion(...
        session_data.metadata.grid_hole, 'SNc');
else
    session_data.metadata.snc_subregion = '';
end

%% Get selected neurons
if isfield(session_data, 'analysis') && isfield(session_data.analysis, 'selected_neurons')
    selected_neurons = session_data.analysis.selected_neurons;
else
    nClusters = height(session_data.spikes.cluster_info);
    selected_neurons = false(nClusters, 1);
    fprintf('Warning: No screening results. All neurons marked NOT SELECTED.\n');
end

%% Create output directory
output_dir = fullfile(project_root, 'figures', 'neuron_summaries');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Generate PDFs for first N neurons
nClusters = height(session_data.spikes.cluster_info);
n_to_process = min(MAX_NEURONS, nClusters);

fprintf('\nGenerating PDFs for %d neurons (of %d total)...\n\n', n_to_process, nClusters);

for i = 1:n_to_process
    fprintf('Processing neuron %d/%d...\n', i, n_to_process);
    tic;
    generate_neuron_summary_pdf(session_data, selected_neurons, ...
        TEST_SESSION, output_dir, 'ClusterIndex', i);
    fprintf('  Completed in %.1f seconds.\n', toc);
end

fprintf('\n=== TEST COMPLETE ===\n');
fprintf('Output directory: %s\n', output_dir);
fprintf('Check the generated PDFs to verify layout.\n');
