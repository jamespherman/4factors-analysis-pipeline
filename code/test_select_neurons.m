%% Test script for the analysis.select_neurons controller.
%
% This script verifies that the select_neurons function correctly dispatches
% to the appropriate area-specific screening function and returns a valid
% list of neuron IDs.
%
% Test Cases:
% 1. SC: Verifies the 'SC' brain area case.
% 2. SNc: Verifies the 'SNc' brain area case.
% 3. Default: Verifies the default case for any other brain area.
%

%% Setup
clear; clc; close all;

% Add the code directory to the MATLAB path to access functions
addpath(fileparts(mfilename('fullpath')));

fprintf('===== Running test_select_neurons =====\n\n');

% Programmatically find data paths
try
    onedrive_path = utils.findOneDrive();
    if isempty(onedrive_path)
        error('OneDrive path not found. Please check utils.findOneDrive.');
    end

    % Define paths for manifest and data directory
    data_dir = fullfile(onedrive_path, 'data', 'processed_data');
    manifest_path = fullfile(data_dir, 'session_manifest.csv');

    % Read the manifest to find test session IDs
    manifest = readtable(manifest_path);

    % Find the first 'SC' session
    sc_session_row = find(strcmp(manifest.brain_area, 'SC'), 1, 'first');
    if isempty(sc_session_row)
        error('No SC session found in the manifest.');
    end
    sc_test_id = manifest.unique_id{sc_session_row};

    % Find the first 'SNc' session
    snc_session_row = find(strcmp(manifest.brain_area, 'SNc'), 1, 'first');
    if isempty(snc_session_row)
        error('No SNc session found in the manifest.');
    end
    snc_test_id = manifest.unique_id{snc_session_row};

catch ME
    fprintf('[FAIL] Setup failed: %s\n', ME.message);
    return; % Stop execution if setup fails
end

fprintf('Setup complete. Found test IDs:\n');
fprintf('  - SC: %s\n', sc_test_id);
fprintf('  - SNc: %s\n\n', snc_test_id);


%% Test 1: SC Selection Case
try
    fprintf('--- Test 1: SC Selection ---\n');

    % Load session data for the SC case
    session_data_sc = data_handling.load_session(sc_test_id, manifest_path, data_dir);

    % Call the controller function to select neurons
    selected_ids_sc = analysis.select_neurons(session_data_sc);

    % --- Verification ---
    % 1. Check if the output is a numeric column vector
    assert(isnumeric(selected_ids_sc) && iscolumn(selected_ids_sc), ...
        'Output is not a numeric column vector.');

    % 2. Check that all selected IDs are part of the original neuron list
    original_ids_sc = session_data_sc.spikes.cluster_info.cluster_id;
    assert(all(ismember(selected_ids_sc, original_ids_sc)), ...
        'Selected IDs are not a subset of original IDs.');

    fprintf('[PASS] SC test ran successfully. Selected %d out of %d total neurons.\n\n', ...
            numel(selected_ids_sc), height(session_data_sc.spikes.cluster_info));

catch ME
    fprintf('[FAIL] SC test failed: %s\n\n', ME.message);
end


%% Test 2: SNc Selection Case
try
    fprintf('--- Test 2: SNc Selection ---\n');

    % Load session data for the SNc case
    session_data_snc = data_handling.load_session(snc_test_id, manifest_path, data_dir);

    % Call the controller function
    selected_ids_snc = analysis.select_neurons(session_data_snc);

    % --- Verification ---
    % 1. Check if the output is a numeric column vector
    assert(isnumeric(selected_ids_snc) && iscolumn(selected_ids_snc), ...
        'Output is not a numeric column vector.');

    % 2. Check that all selected IDs are part of the original neuron list
    original_ids_snc = session_data_snc.spikes.cluster_info.cluster_id;
    assert(all(ismember(selected_ids_snc, original_ids_snc)), ...
        'Selected IDs are not a subset of original IDs.');

    fprintf('[PASS] SNc test ran successfully. Selected %d out of %d total neurons.\n\n', ...
            numel(selected_ids_snc), height(session_data_snc.spikes.cluster_info));

catch ME
    fprintf('[FAIL] SNc test failed: %s\n\n', ME.message);
end


%% Test 3: Default/Other Area Case
try
    fprintf('--- Test 3: Default (Other Area) Selection ---\n');

    % Use the loaded SC data as a base and modify the brain area
    session_data_other = session_data_sc;
    session_data_other.metadata.brain_area = 'Cortex'; % Unrecognized area

    % Call the controller, which should trigger the default case
    selected_ids_other = analysis.select_neurons(session_data_other);

    % --- Verification ---
    % The function should return all neurons marked as 'good'
    expected_ids = utils.get_good_neurons(session_data_other);

    % 1. Check if the selected IDs are identical to the expected 'good' IDs
    % Use sort to ensure the order doesn't matter for comparison
    assert(isequal(sort(selected_ids_other), sort(expected_ids)), ...
        'Default case did not return the list of "good" neurons.');

    fprintf('[PASS] Default case handled correctly. Selected %d "good" neurons.\n\n', ...
            numel(selected_ids_other));

catch ME
    fprintf('[FAIL] Default case test failed: %s\n\n', ME.message);
end

fprintf('===== Test script finished =====\n');
