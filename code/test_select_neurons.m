%% Test script for the analysis.select_neurons function
%
% This script verifies that the analysis.select_neurons controller
% correctly calls the appropriate area-specific screening function and
% returns a valid list of neuron IDs. It tests the logic for SC, SNc,
% and an unrecognized brain area.

%% Standard Setup
% This block assumes the script is run from the project's root directory.
% It defines key path variables and adds necessary folders to the MATLAB path.

% Define project root and data directory paths.
project_root = pwd;
data_base_path = utils.findOneDrive(); % Standard way to find the data directory.

% Add all sub-directories within the 'code' folder to the MATLAB path.
addpath(genpath(fullfile(project_root, 'code')));

% Construct the full path to the session manifest.
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');


%% Test-Specific Setup
% Read the session manifest to find session IDs for testing.
manifest = readtable(manifest_path);

% Find a unique_id for an 'SC' session.
scRows = find(strcmp(manifest.brain_area, 'SC'));
if isempty(scRows)
    error('No SC sessions found in the manifest.');
end
scTestId = manifest.unique_id{scRows(1)};

% Find a unique_id for an 'SNc' session.
sncRows = find(strcmp(manifest.brain_area, 'SNc'));
if isempty(sncRows)
    error('No SNc sessions found in the manifest.');
end
sncTestId = manifest.unique_id{sncRows(1)};

fprintf('Setup complete. Using ID %s for SC and ID %s for SNc.\n', ...
    scTestId, sncTestId);


%% Test 1: SC Selection Case
try
    % Load the SC session data using the standardized function.
    sessionDataSc = data_handling.load_session(scTestId, manifest_path, data_base_path);

    % Call the neuron selection function.
    selectedIdsSc = analysis.select_neurons(sessionDataSc);

    % --- Verification ---
    % 1. Check that the output is a numeric column vector.
    assert(isnumeric(selectedIdsSc) && iscolumn(selectedIdsSc), ...
        'Output is not a numeric column vector.');

    % 2. Check that all selected IDs are part of the original neuron list.
    originalIdsSc = sessionDataSc.spikes.cluster_info.cluster_id;
    assert(all(ismember(selectedIdsSc, originalIdsSc)), ...
        'Selected IDs are not a subset of original IDs.');

    % --- Report Results ---
    fprintf('[PASS] SC test ran successfully. Selected %d out of %d total neurons.\n', ...
        numel(selectedIdsSc), height(sessionDataSc.spikes.cluster_info));
catch ME
    fprintf('[FAIL] SC test failed: %s\n', ME.message);
end


%% Test 2: SNc Selection Case
try
    % Load the SNc session data using the standardized function.
    sessionDataSnc = data_handling.load_session(sncTestId, manifest_path, data_base_path);

    % Call the neuron selection function.
    selectedIdsSnc = analysis.select_neurons(sessionDataSnc);

    % --- Verification ---
    % 1. Check that the output is a numeric column vector.
    assert(isnumeric(selectedIdsSnc) && iscolumn(selectedIdsSnc), ...
        'Output is not a numeric column vector.');

    % 2. Check that all selected IDs are part of the original neuron list.
    originalIdsSnc = sessionDataSnc.spikes.cluster_info.cluster_id;
    assert(all(ismember(selectedIdsSnc, originalIdsSnc)), ...
        'Selected IDs are not a subset of original IDs.');

    % --- Report Results ---
    fprintf('[PASS] SNc test ran successfully. Selected %d out of %d total neurons.\n', ...
        numel(selectedIdsSnc), height(sessionDataSnc.spikes.cluster_info));
catch ME
    fprintf('[FAIL] SNc test failed: %s\n', ME.message);
end


%% Test 3: Default/Other Area Case
try
    % Use the loaded SC data as a base and modify the brain area.
    % This test does not require loading new data, so it's fine as is.
    sessionDataOther = sessionDataSc;
    sessionDataOther.metadata.brain_area = 'Cortex'; % An unrecognized area.

    % Call the neuron selection function.
    selectedIdsOther = analysis.select_neurons(sessionDataOther);

    % --- Verification ---
    % The function should return an empty list for unrecognized areas.
    assert(isempty(selectedIdsOther), ...
        'The function did not return an empty list for an unrecognized area.');

    % --- Report Results ---
    fprintf('[PASS] Default case for unrecognized brain area handled correctly.\n');
catch ME
    fprintf('[FAIL] Default case test failed: %s\n', ME.message);
end
