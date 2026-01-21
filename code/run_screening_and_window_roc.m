%% run_screening_and_window_roc.m
%
% Focused pipeline script that runs ONLY:
%   1. Core data preparation (for ALL neurons)
%   2. Refined neuron screening (FR threshold + task modulation)
%   3. Window-based ROC analysis
%   4. Per-neuron diagnostic PDFs (for ALL neurons)
%
% This script skips: baseline comparison, bin-by-bin ROC, ANOVA,
% behavioral analyses, and population decoding.
%
% Author: Claude Code
% Date: 2026-01-18

%% Setup
clear; clc; close all;

% --- USER TOGGLES ---
force_rerun = struct(...
    'screening', false, ...   % Set to true to recompute screening
    'window_roc', true, ...  % Set to true to recompute window ROC
    'diag_pdfs', false ...    % Set to true to regenerate PDFs
);

% Optional: specify specific sessions to process (empty = all sessions)
sessions_to_process = {};  % Empty = process all sessions
% --- END USER TOGGLES ---

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

tic;
giveFeed = @(x) disp([num2str(round(toc, 1)) 's - ' x]);

%% Load Manifest
giveFeed('Loading session manifest...');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
if ~exist(manifest_path, 'file')
    error('Session manifest not found at: %s', manifest_path);
end
manifest = readtable(manifest_path);
giveFeed(sprintf('Manifest loaded: %d sessions.', height(manifest)));

%% Load Analysis Plan
giveFeed('Loading analysis plan...');
[~, analysis_plan] = define_task_conditions();
giveFeed('Analysis plan loaded.');

%% Filter sessions if specified
if ~isempty(sessions_to_process)
    session_mask = ismember(manifest.unique_id, sessions_to_process);
    if sum(session_mask) == 0
        error('None of the specified sessions found in manifest.');
    end
    manifest = manifest(session_mask, :);
    giveFeed(sprintf('Filtered to %d specified sessions.', height(manifest)));
end

%% Iterate Through Sessions
for i = 1:height(manifest)
    session_id = manifest.unique_id{i};
    giveFeed(sprintf('\n=== Processing session %d/%d: %s ===', i, height(manifest), session_id));

    %% Load Session Data
    local_processed_dir = fullfile(project_root, 'data', ...
        'processed', session_id);
    local_session_data_path = fullfile(local_processed_dir, ...
        [session_id, '_session_data.mat']);

    if exist(local_session_data_path, 'file')
        giveFeed('Loading cached data...');
        load(local_session_data_path, 'session_data');
    else
        % Load from OneDrive
        one_drive_path = findOneDrive;
        source_path = fullfile(one_drive_path, 'Neuronal Data Analysis', ...
            session_id, [session_id '_session_data.mat']);

        if ~exist(source_path, 'file')
            warning('Session data not found for %s. Skipping.', session_id);
            continue;
        end

        giveFeed('Loading from OneDrive...');
        load(source_path, 'session_data');
    end

    % Add metadata from manifest
    session_data.metadata = table2struct(manifest(i, :));
    session_data.metadata.unique_id = session_id;

    % Parse SNc subregion
    session_data.metadata.snc_subregion = parse_snc_subregion(...
        session_data.metadata.grid_hole, session_data.metadata.brain_area);

    data_updated = false;

    %% On-Demand Metrics Calculation
    if ~isfield(session_data, 'metrics') || ~isfield(session_data.metrics, 'baseline_frs')
        giveFeed('Computing metrics on-the-fly...');
        codes = initCodes();
        config = analysis_plan.baseline_fr_config;
        task_code = codes.(['uniqueTaskCode_' config.task_name]);
        trial_mask = session_data.trialInfo.taskCode == task_code;
        baseline_fr = calculate_baseline_fr(session_data, config.event, config.window, trial_mask);
        session_data.metrics.baseline_frs = baseline_fr;

        % Waveform metrics
        nClusters = numel(session_data.spikes.cluster_info.cluster_id);
        for i_cluster = 1:nClusters
            mean_waveform = session_data.spikes.wfMeans{i_cluster};
            [~, max_var_chan] = max(var(mean_waveform, [], 2));
            session_data.metrics.wf_metrics(i_cluster, 1) = ...
                calculate_waveform_metrics(mean_waveform(max_var_chan, :), 30000);
        end
        data_updated = true;
        giveFeed('Metrics computed.');
    end

    %% Define conditions
    conditions = define_task_conditions(session_data);

    %% Step 1: Core Data Preparation (ALL neurons)
    nClusters = numel(session_data.spikes.cluster_info.cluster_id);
    all_neurons_mask = true(nClusters, 1);
    alignment_events = analysis_plan.events;

    needs_core_data = force_rerun.screening || force_rerun.window_roc || ...
        ~isfield(session_data, 'analysis') || ...
        ~isfield(session_data.analysis, 'core_data');

    if needs_core_data
        giveFeed(sprintf('Preparing core_data for %d neurons...', nClusters));
        core_data = prepare_core_data(session_data, all_neurons_mask, ...
            alignment_events, analysis_plan);
        session_data.analysis.core_data = core_data;
        data_updated = true;
        giveFeed('Core data prepared.');
    else
        giveFeed('Loading existing core_data...');
        core_data = session_data.analysis.core_data;
    end

    %% Step 2: Neuron Screening
    needs_screening = force_rerun.screening || ...
        ~isfield(session_data.analysis, 'screening_info');

    if needs_screening
        giveFeed('Running refined neuron screening...');
        [selected_neurons, screening_info] = apply_neuron_screening(...
            session_data, core_data, analysis_plan);

        session_data.analysis.selected_neurons = selected_neurons;
        session_data.analysis.screening_info = screening_info;

        % Determine scSide for SC sessions
        if strcmp(session_data.metadata.brain_area, 'SC')
            grid_hole_str = session_data.metadata.grid_hole;
            coords = sscanf(grid_hole_str, '(%f, %f)');
            if ~isempty(coords)
                if coords(1) < 0
                    scSide = 'left';
                elseif coords(1) > 0
                    scSide = 'right';
                else
                    scSide = 'unknown';
                end
            else
                scSide = 'unknown';
            end
            session_data.analysis.scSide = scSide;
            giveFeed(sprintf('SC Side: %s', scSide));
        end

        data_updated = true;

        % Print screening summary
        n_included = sum(selected_neurons);
        giveFeed(sprintf('Screening: %d/%d neurons included (%.1f%%)', ...
            n_included, nClusters, 100 * n_included / nClusters));
    else
        giveFeed('Loading existing screening results...');
        selected_neurons = session_data.analysis.selected_neurons;
        screening_info = session_data.analysis.screening_info;
    end

    %% Step 3: Window-Based ROC Analysis
    needs_window_roc = force_rerun.window_roc || ...
        ~isfield(session_data.analysis, 'window_roc');

    if needs_window_roc && isfield(analysis_plan, 'window_roc_plan')
        giveFeed('Running window-based ROC analysis...');

        % Compute preferred locations for SC neurons (if using 'preferred' mode)
        preferred_locations = [];
        if strcmp(session_data.metadata.brain_area, 'SC') && ...
           isfield(analysis_plan.window_roc_plan, 'location_mode') && ...
           strcmp(analysis_plan.window_roc_plan.location_mode, 'preferred')
            giveFeed('Computing per-neuron preferred locations...');
            preferred_locations = compute_preferred_locations(...
                session_data, core_data, conditions, analysis_plan);
            session_data.analysis.preferred_locations = preferred_locations;

            % Report distribution
            pref_locs = preferred_locations.for_factors;
            for loc = 1:4
                n_at_loc = sum(pref_locs == loc);
                giveFeed(sprintf('  Location %d: %d neurons', loc, n_at_loc));
            end
        end

        % Run window ROC analysis
        window_roc_results = analyze_window_roc(session_data, ...
            conditions, analysis_plan, core_data, preferred_locations);
        session_data.analysis.window_roc = window_roc_results;

        data_updated = true;
        giveFeed('Window ROC analysis complete.');
    else
        giveFeed('Window ROC already complete or skipped.');
    end

    %% Step 4: Generate PDFs for ALL neurons
    diag_output_dir = fullfile(project_root, 'figures', session_id);

    needs_pdfs = force_rerun.diag_pdfs || force_rerun.screening || ...
        ~exist(diag_output_dir, 'dir') || ...
        isempty(dir(fullfile(diag_output_dir, '*.pdf')));

    if needs_pdfs
        giveFeed(sprintf('Generating PDFs for %d neurons...', nClusters));
        if ~exist(diag_output_dir, 'dir')
            mkdir(diag_output_dir);
        end

        generate_neuron_summary_pdf(session_data, selected_neurons, ...
            session_id, diag_output_dir, 'ScreeningInfo', screening_info);
        giveFeed('PDF generation complete.');
    else
        giveFeed('PDFs already exist, skipping...');
    end

    %% Save Updated Data
    if data_updated
        giveFeed('Saving updated session data...');
        if ~exist(local_processed_dir, 'dir')
            mkdir(local_processed_dir);
        end
        save(local_session_data_path, 'session_data', '-v7.3');
        giveFeed('Save complete.');
    end

    giveFeed(sprintf('=== Finished %s ===\n', session_id));
end

%% Summary
giveFeed('');
giveFeed('========================================');
giveFeed('ALL SESSIONS COMPLETE');
giveFeed('========================================');
giveFeed(sprintf('Total time: %.1f seconds', toc));
giveFeed('');
giveFeed('Next steps:');
giveFeed('1. Check figures/{session_id}/ for neuron PDFs');
giveFeed('2. Run aggregate_analysis_results.m to generate summary figures');
