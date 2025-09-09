%% Test script for the analysis.select_neurons function
%
% This script verifies that the analysis.select_neurons controller
% correctly calls the appropriate area-specific screening function and
% returns a valid list of neuron IDs. It tests the logic for SC, SNc,
% and an unrecognized brain area.

%% Setup
% Add the code directory to the MATLAB path to ensure all functions are
% accessible.
addpath(fileparts(mfilename('fullpath')));

% Find the path to the manifest and the data directory.
[~, datapath] = utils.findOneDrive();
manifestPath = fullfile(datapath, '..', 'config', 'session_manifest.csv');

% Read the session manifest.
manifest = readtable(manifestPath);

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
    % Load the SC session data.
    sessionDataSc = data_handling.load_session(scTestId, datapath);

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
    fprintf(['[PASS] SC test ran successfully. Selected %d out of %d ' ...
        'total neurons.\n'], numel(selectedIdsSc), ...
        height(sessionDataSc.spikes.cluster_info));
catch ME
    fprintf('[FAIL] SC test failed: %s\n', ME.message);
end


%% Test 2: SNc Selection Case
try
    % Load the SNc session data.
    sessionDataSnc = data_handling.load_session(sncTestId, datapath);

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
    fprintf(['[PASS] SNc test ran successfully. Selected %d out of %d ' ...
        'total neurons.\n'], numel(selectedIdsSnc), ...
        height(sessionDataSnc.spikes.cluster_info));
catch ME
    fprintf('[FAIL] SNc test failed: %s\n', ME.message);
end


%% Test 3: Default/Other Area Case
try
    % Use the loaded SC data as a base and modify the brain area.
    sessionDataOther = sessionDataSc;
    sessionDataOther.metadata.brain_area = 'Cortex'; % An unrecognized area.

    % Call the neuron selection function.
    selectedIdsOther = analysis.select_neurons(sessionDataOther);

    % --- Verification ---
    % The function should return an empty list for unrecognized areas.
    assert(isempty(selectedIdsOther), ...
        'The function did not return an empty list for an unrecognized area.');

    % --- Report Results ---
    fprintf(['[PASS] Default case for unrecognized brain area handled ' ...
        'correctly.\n']);
catch ME
    fprintf('[FAIL] Default case test failed: %s\n', ME.message);
end
