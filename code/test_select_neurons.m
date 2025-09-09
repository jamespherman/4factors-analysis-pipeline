%% Test script for the analysis.select_neurons function
%
% This script verifies that the analysis.select_neurons controller
% correctly calls the appropriate area-specific screening function and
% returns a valid list of neuron IDs. It tests the logic for SC, SNc,
% and an unrecognized brain area.

%% Setup
% Add the code directory to the MATLAB path to ensure all functions are accessible.
addpath(fileparts(mfilename('fullpath')));

% Find the path to the manifest and the data directory.
[~, datapath] = utils.findOneDrive();
manifest_path = fullfile(datapath, '..', 'config', 'session_manifest.csv');

% Read the session manifest.
manifest = readtable(manifest_path);

% Find a unique_id for an 'SC' session.
sc_rows = find(strcmp(manifest.brain_area, 'SC'));
if isempty(sc_rows)
    error('No SC sessions found in the manifest.');
end
sc_test_id = manifest.unique_id{sc_rows(1)};

% Find a unique_id for an 'SNc' session.
snc_rows = find(strcmp(manifest.brain_area, 'SNc'));
if isempty(snc_rows)
    error('No SNc sessions found in the manifest.');
end
snc_test_id = manifest.unique_id{snc_rows(1)};

fprintf('Setup complete. Using ID %s for SC and ID %s for SNc.\n', sc_test_id, snc_test_id);


%% Test 1: SC Selection Case
try
    % Load the SC session data.
    session_data_sc = data_handling.load_session(sc_test_id, datapath);

    % Call the neuron selection function.
    selected_ids_sc = analysis.select_neurons(session_data_sc);

    % --- Verification ---
    % 1. Check that the output is a numeric column vector.
    assert(isnumeric(selected_ids_sc) && iscolumn(selected_ids_sc), 'Output is not a numeric column vector.');

    % 2. Check that all selected IDs are part of the original neuron list.
    original_ids_sc = session_data_sc.spikes.cluster_info.cluster_id;
    assert(all(ismember(selected_ids_sc, original_ids_sc)), 'Selected IDs are not a subset of original IDs.');

    % --- Report Results ---
    fprintf('[PASS] SC test ran successfully. Selected %d out of %d total neurons.\n', ...
            numel(selected_ids_sc), height(session_data_sc.spikes.cluster_info));
catch ME
    fprintf('[FAIL] SC test failed: %s\n', ME.message);
end


%% Test 2: SNc Selection Case
try
    % Load the SNc session data.
    session_data_snc = data_handling.load_session(snc_test_id, datapath);

    % Call the neuron selection function.
    selected_ids_snc = analysis.select_neurons(session_data_snc);

    % --- Verification ---
    % 1. Check that the output is a numeric column vector.
    assert(isnumeric(selected_ids_snc) && iscolumn(selected_ids_snc), 'Output is not a numeric column vector.');

    % 2. Check that all selected IDs are part of the original neuron list.
    original_ids_snc = session_data_snc.spikes.cluster_info.cluster_id;
    assert(all(ismember(selected_ids_snc, original_ids_snc)), 'Selected IDs are not a subset of original IDs.');

    % --- Report Results ---
    fprintf('[PASS] SNc test ran successfully. Selected %d out of %d total neurons.\n', ...
            numel(selected_ids_snc), height(session_data_snc.spikes.cluster_info));
catch ME
    fprintf('[FAIL] SNc test failed: %s\n', ME.message);
end


%% Test 3: Default/Other Area Case
try
    % Use the loaded SC data as a base and modify the brain area.
    session_data_other = session_data_sc;
    session_data_other.metadata.brain_area = 'Cortex'; % An unrecognized area.

    % Call the neuron selection function.
    selected_ids_other = analysis.select_neurons(session_data_other);

    % --- Verification ---
    % The function should return an empty list for unrecognized areas.
    assert(isempty(selected_ids_other), 'The function did not return an empty list for an unrecognized brain area.');

    % --- Report Results ---
    fprintf('[PASS] Default case for unrecognized brain area handled correctly.\n');
catch ME
    fprintf('[FAIL] Default case test failed: %s\n', ME.message);
end
