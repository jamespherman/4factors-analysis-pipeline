% run_4factors_analysis.m
%
% This script serves as the top-level controller for the 4-factors task
% analysis pipeline. It processes multiple sessions listed in a manifest
% file, handling neuron selection and data preparation in an idempotent
% manner.
%
% The script performs the following main steps for each session:
%   1. On-Demand Neuron Selection: Selects neurons if not already done.
%   2. On-Demand Data Preparation: Prepares binned spike data for key
%      events if not already prepared.
%   3. Run Analyses: (Placeholder for future implementation).
%
% This script is designed to be run from the root of the project directory.
%
% See also:
%   data_handling.load_session, analysis.select_neurons,
%   utils.alignAndBinSpikes, config/session_manifest.csv

% --- Script Initialization ---
% Clear workspace and close all figures to ensure a clean run.
clear;
close all;

fprintf('Starting the 4-factors analysis pipeline...\n');

%% --- 1. Standard Setup ---
% This section defines the necessary paths for the analysis pipeline. It
% assumes the script is run from the root of the project directory.

fprintf('1. Setting up paths...\n');

% Get the project root directory
projectRoot = pwd;

% Add the 'code' directory to the MATLAB path
codeDir = fullfile(projectRoot, 'code');
addpath(codeDir);

% Define the path to the session manifest file
manifestPath = fullfile(projectRoot, 'config', 'session_manifest.csv');

% Verify that the manifest file exists
if ~isfile(manifestPath)
    error('Session manifest file not found at: %s', manifestPath);
end

%% --- 2. Load Manifest and Loop Through Sessions ---
fprintf('2. Loading session manifest and finding data directory...\n');

% Load the session manifest
manifest = readtable(manifestPath, 'PreserveVariableNames', true);

% Get the base path for the data directory on OneDrive
try
    dataBasePath = utils.findOneDrive();
catch ME
    error('Could not find OneDrive base path. Error: %s', ME.message);
end

fprintf('Found data directory at: %s\n', dataBasePath);

% Loop through each session in the manifest
for i = 1:height(manifest)
    sessionInfo = manifest(i, :);
    uniqueId = sessionInfo.unique_id{1};

    fprintf('\n--- Processing session: %s ---\n', uniqueId);

    try
        % Load the session data once at the beginning of the loop
        fprintf('   -> Loading session_data...\n');
        session_data = data_handling.load_session(uniqueId, manifestPath, dataBasePath);
        dataNeedsSaving = false; % Flag to track if session_data was modified

        %% --- 3. Step 1: On-Demand Neuron Selection (Idempotent) ---
        if strcmpi(sessionInfo.selection_status, 'pending')
            fprintf('Status is "pending". Running neuron selection...\n');

            selectedIds = analysis.select_neurons(session_data);
            session_data.analysis.selected_neurons = selectedIds;
            fprintf('   -> Found %d selected neurons.\n', numel(selectedIds));
            dataNeedsSaving = true;

            % Update status in the manifest and save immediately for robustness
            manifest.selection_status{i} = 'complete';
            writetable(manifest, manifestPath);
            fprintf('   -> Status updated to "complete" and manifest saved.\n');
        else
            fprintf('Status is "%s". Skipping neuron selection.\n', ...
                sessionInfo.selection_status{1});
        end

        %% --- 4. Step 2: On-Demand Data Preparation (Idempotent) ---
        % Re-read the session info in case the status was just updated
        sessionInfo = manifest(i, :);
        if strcmpi(sessionInfo.dataprep_status, 'pending')
            fprintf('Status is "pending". Running data preparation...\n');

            if ~isfield(session_data, 'analysis') || ~isfield(session_data.analysis, 'selected_neurons')
                fprintf('ERROR: Neuron selection has not been run for %s. Skipping data prep.\n', uniqueId);
                continue;
            end

            % Get analysis plan from the dedicated function
            analysisPlan = analysis.define_analysis_plan();
            alignEvents = analysisPlan.alignEvents;
            timeWindow = analysisPlan.timeWindow;
            binWidth = analysisPlan.binWidth;

            fprintf('   -> Preparing data for events: %s\n', strjoin(alignEvents, ', '));

            preparedData = struct();
            selectedNeuronIds = session_data.analysis.selected_neurons;
            nNeurons = numel(selectedNeuronIds);

            for evIdx = 1:numel(alignEvents)
                eventName = alignEvents{evIdx};
                fprintf('      - Aligning to event: %s\n', eventName);

                eventTimes = session_data.trialInfo.timing.(eventName);
                nTrials = numel(eventTimes);
                nBins = floor((timeWindow(2) - timeWindow(1)) / binWidth);
                binnedSpikesEvent = zeros(nNeurons, nTrials, nBins);

                for nIdx = 1:nNeurons
                    neuronId = selectedNeuronIds(nIdx);
                    isNeuronSpike = (session_data.spikes.clusters == neuronId);
                    neuronSpikeTimes = session_data.spikes.times(isNeuronSpike);

                    [~, binnedCounts] = utils.alignAndBinSpikes(neuronSpikeTimes, ...
                        eventTimes, timeWindow(1), timeWindow(2), binWidth);

                    binnedSpikesEvent(nIdx, :, :) = binnedCounts;
                end
                preparedData.(eventName) = binnedSpikesEvent;
            end

            session_data.analysis.prepared_data = preparedData;
            dataNeedsSaving = true;

            % Update status in the manifest and save immediately
            manifest.dataprep_status{i} = 'complete';
            writetable(manifest, manifestPath);
            fprintf('   -> Status updated to "complete" and manifest saved.\n');
        else
            fprintf('Status is "%s". Skipping data preparation.\n', ...
                sessionInfo.dataprep_status{1});
        end

        %% --- Save session_data if it was modified ---
        if dataNeedsSaving
            fprintf('   -> Saving updated session_data to disk...\n');
            save(fullfile(dataBasePath, uniqueId, [uniqueId, '_session_data.mat']), 'session_data', '-v7.3');
            fprintf('   -> Save complete.\n');
        end

        %% --- 5. Step 3: Run Analyses (Future Placeholder) ---
        fprintf('--- Session %s processing complete. ---\n', uniqueId);

    catch ME
        fprintf('!!! ERROR processing session %s: %s\n', uniqueId, ME.message);
        fprintf('Continuing to next session.\n');
    end
end

fprintf('\nAnalysis pipeline finished.\n');
