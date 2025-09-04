%% Test script for data_handling.load_session

% This script is designed to be run from the project's root directory.
fprintf('Starting test for data_handling.load_session.m\n');

% Add the 'code' directory to the MATLAB path
try
    addpath(genpath('code'));
    fprintf('[INFO] Added "code" directory to MATLAB path.\n');
catch ME
    fprintf('[ERROR] Could not add "code" directory to path. %s\n', ME.message);
    return;
end

%% --- Test 1: Successful Data Loading ---
fprintf('\n--- Running Test 1: Successful Data Loading ---\n');
try
    % Programmatically determine paths
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    project_root = fileparts(script_dir);
    manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
    fprintf('[INFO] Manifest path: %s\n', manifest_path);

    % Get OneDrive path and construct data_base_path
    oneDrivePath = findOneDrive();
    data_base_path = fullfile(oneDrivePath, 'Neuronal Data Analysis');
    fprintf('[INFO] Data base path: %s\n', data_base_path);

    % Read the manifest to get a valid unique_id
    manifest_table = readtable(manifest_path);
    TEST_UNIQUE_ID = manifest_table.unique_id{1};
    fprintf('[INFO] Using test unique_id: %s\n', TEST_UNIQUE_ID);

    % Call the function to load session data
    session_data = data_handling.load_session(TEST_UNIQUE_ID, manifest_path, data_base_path);

    % Verification
    assert(isfield(session_data, 'metadata'));
    fprintf('[PASS] session_data contains the ".metadata" field.\n');

    assert(isfield(session_data.metadata, 'unique_id'));
    fprintf('[PASS] session_data.metadata contains the "unique_id" field.\n');

    assert(strcmp(session_data.metadata.unique_id, TEST_UNIQUE_ID));
    fprintf('[PASS] Metadata''s unique_id matches the test ID.\n');

    % A simple check for some other expected field (e.g., 'spikes')
    % This assumes 'spikes' is a field in the original data.
    % If not, this check might need to be adapted based on actual data structure.
    assert(isfield(session_data, 'spikes'));
    fprintf('[PASS] Original data fields (e.g., "spikes") are intact.\n');

    fprintf('--- Test 1 COMPLETE ---\n');

catch ME
    fprintf('[FAIL] Test 1 threw an unexpected error: %s\n', ME.message);
    fprintf('--- Test 1 FAILED ---\n');
end


%% --- Test 2: Graceful Failure for Non-Existent Session ---
fprintf('\n--- Running Test 2: Graceful Failure for Non-Existent Session ---\n');
FAKE_UNIQUE_ID = 'This_ID_Does_Not_Exist_123';
try
    % This call should fail and throw an error
    data_handling.load_session(FAKE_UNIQUE_ID, manifest_path, data_base_path);

    % If it reaches here, the error was not thrown, which is a failure
    fprintf('[FAIL] The function did not throw an error for a fake unique_id.\n');

catch ME
    % The error is expected, so this is a pass
    fprintf('[PASS] The function correctly threw an error for a non-existent unique_id: %s\n', ME.message);
end
fprintf('--- Test 2 COMPLETE ---\n');


fprintf('\nAll tests finished.\n');
