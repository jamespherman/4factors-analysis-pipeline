%% generate_all_neuron_pdfs.m
%
% Standalone batch script to generate neuron summary PDFs for all neurons
% across all completed sessions. This script does NOT require the full
% analysis pipeline to have run first - it only needs the session_data.mat
% files to exist.
%
% The script:
%   1. Reads session_manifest.csv
%   2. Loops through all sessions with analysis_status == 'complete'
%   3. Loads each _session_data.mat file
%   4. Generates PDFs for ALL neurons (not just selected ones)
%   5. Saves to figures/neuron_summaries/
%   6. Reports progress and handles errors gracefully
%
% Output naming: {unique_id}_{cluster_id:03d}_summary.pdf
%
% Author: Claude (with James Herman)
% Date: 2025-01-08

%% Setup
clear; clc; close all;

% Setup paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

% Timing
tic;
giveFeed = @(x) fprintf('[%6.1fs] %s\n', toc, x);

%% Configuration
% Set to true to regenerate PDFs even if they already exist
FORCE_REGENERATE = false;

% Set to true to skip sessions that encounter errors (otherwise stops)
CONTINUE_ON_ERROR = true;

%% Load Manifest
giveFeed('Loading session manifest...');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
if ~exist(manifest_path, 'file')
    error('generate_all_neuron_pdfs:manifestNotFound', ...
        'Session manifest not found at: %s', manifest_path);
end
manifest = readtable(manifest_path);
giveFeed(sprintf('Manifest loaded: %d sessions total.', height(manifest)));

%% Filter to completed sessions
% Only process sessions where analysis_status == 'complete'
completed_mask = strcmp(manifest.analysis_status, 'complete');
sessions_to_process = manifest(completed_mask, :);
n_sessions = height(sessions_to_process);

if n_sessions == 0
    warning('No sessions with analysis_status == ''complete'' found.');
    return;
end

giveFeed(sprintf('Found %d completed sessions to process.', n_sessions));

%% Create output directory
output_dir = fullfile(project_root, 'figures', 'neuron_summaries');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    giveFeed(sprintf('Created output directory: %s', output_dir));
end

%% Initialize counters
total_neurons_processed = 0;
total_neurons_skipped = 0;
sessions_with_errors = {};

%% Process each session
for i_session = 1:n_sessions
    session_id = sessions_to_process.unique_id{i_session};

    giveFeed(sprintf('=== Processing session %d/%d: %s ===', ...
        i_session, n_sessions, session_id));

    try
        %% Load session data
        % Try local processed directory first
        local_processed_dir = fullfile(project_root, 'data', 'processed', ...
            session_id);
        local_session_data_path = fullfile(local_processed_dir, ...
            [session_id, '_session_data.mat']);

        if exist(local_session_data_path, 'file')
            giveFeed(sprintf('  Loading from local cache...'));
            load(local_session_data_path, 'session_data');
        else
            % Try OneDrive location
            one_drive_path = findOneDrive();
            source_session_data_path = fullfile(one_drive_path, ...
                'Neuronal Data Analysis', session_id, ...
                [session_id '_session_data.mat']);

            if ~exist(source_session_data_path, 'file')
                warning('  Session data not found for %s. Skipping.', session_id);
                sessions_with_errors{end+1} = session_id; %#ok<AGROW>
                continue;
            end

            giveFeed(sprintf('  Loading from OneDrive...'));
            load(source_session_data_path, 'session_data');
        end

        %% Ensure metadata is populated
        session_data.metadata = table2struct(sessions_to_process(i_session, :));

        % Add SNc subregion if applicable
        if strcmp(session_data.metadata.brain_area, 'SNc')
            session_data.metadata.snc_subregion = parse_snc_subregion(...
                session_data.metadata.grid_hole, 'SNc');
        else
            session_data.metadata.snc_subregion = '';
        end

        %% Get neuron count
        nClusters = numel(session_data.spikes.cluster_info.cluster_id);
        giveFeed(sprintf('  Found %d neurons.', nClusters));

        %% Get selected_neurons (if available from screening)
        if isfield(session_data, 'analysis') && ...
                isfield(session_data.analysis, 'selected_neurons')
            selected_neurons = session_data.analysis.selected_neurons;
        else
            % If no screening results, mark all as not selected
            selected_neurons = false(nClusters, 1);
            giveFeed('  Warning: No screening results found. All neurons marked as NOT SELECTED.');
        end

        %% Generate PDFs for each neuron
        for i_neuron = 1:nClusters
            cluster_id = session_data.spikes.cluster_info.cluster_id(i_neuron);

            % Check if PDF already exists
            pdf_filename = sprintf('%s_%03d_summary.pdf', session_id, cluster_id);
            pdf_path = fullfile(output_dir, pdf_filename);

            if exist(pdf_path, 'file') && ~FORCE_REGENERATE
                total_neurons_skipped = total_neurons_skipped + 1;
                continue;
            end

            % Progress report
            fprintf('  Neuron %d/%d (Cluster %d)...', i_neuron, nClusters, cluster_id);

            % Generate PDF for this single neuron
            generate_neuron_summary_pdf(session_data, selected_neurons, ...
                session_id, output_dir, 'ClusterIndex', i_neuron);

            total_neurons_processed = total_neurons_processed + 1;
            fprintf(' done.\n');
        end

        giveFeed(sprintf('  Session complete: %d neurons processed.', nClusters));

    catch ME
        warning('Error processing session %s: %s', session_id, ME.message);
        sessions_with_errors{end+1} = session_id; %#ok<AGROW>

        if ~CONTINUE_ON_ERROR
            rethrow(ME);
        end
    end
end

%% Summary
giveFeed('=== BATCH PROCESSING COMPLETE ===');
fprintf('\n');
fprintf('Sessions processed: %d\n', n_sessions);
fprintf('Neurons with new PDFs: %d\n', total_neurons_processed);
fprintf('Neurons skipped (PDF exists): %d\n', total_neurons_skipped);
fprintf('Sessions with errors: %d\n', length(sessions_with_errors));

if ~isempty(sessions_with_errors)
    fprintf('\nSessions with errors:\n');
    for i = 1:length(sessions_with_errors)
        fprintf('  - %s\n', sessions_with_errors{i});
    end
end

fprintf('\nOutput directory: %s\n', output_dir);
fprintf('Total time: %.1f seconds\n', toc);
