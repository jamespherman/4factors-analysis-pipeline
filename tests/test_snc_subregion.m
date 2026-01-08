%% test_snc_subregion.m
%
% Minimal test script to verify the SNc subregion labeling implementation
% without running the full analysis pipeline.
%
% Tests:
%   1. parse_snc_subregion function with known manifest values
%   2. Condition mask creation via define_task_conditions
%   3. Summary of all manifest sessions with parsed subregions
%
% Usage:
%   >> cd('path/to/4factors-analysis-pipeline')
%   >> run('tests/test_snc_subregion.m')
%
% Author: Claude Code
% Date: 2026-01-08

%% Setup
clear; clc;
fprintf('=== SNc Subregion Labeling Test Suite ===\n\n');

% Add utils to path
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(fullfile(project_root, 'code', 'utils'));

test_passed = 0;
test_failed = 0;

%% Test 1: parse_snc_subregion with known values
fprintf('--- Test 1: parse_snc_subregion Direct Tests ---\n');

% Test case structure: {grid_hole, brain_area, expected_result, description}
test_cases = {
    '(+5, +3)', 'SNc', 'rvmSNc', 'Y=3 (boundary) should be rvmSNc';
    '(+5, +1)', 'SNc', 'cdlSNc', 'Y=1 should be cdlSNc';
    '(+4, +2)', 'SC',  '',       'SC sessions should have empty subregion';
    '(+3, +4)', 'SNc', 'rvmSNc', 'Y=4 should be rvmSNc';
    '(+2, +2)', 'SNc', 'cdlSNc', 'Y=2 should be cdlSNc';
    '(+6, +2)', 'SNc', 'cdlSNc', 'Y=2 should be cdlSNc';
    '(-2,+1.5)','SNc', 'cdlSNc', 'Y=1.5 (no space) should be cdlSNc';
    '(+3.5, +2.5)', 'SNc', 'cdlSNc', 'Y=2.5 (decimal) should be cdlSNc';
    '(+4,+0.5)', 'SC', '', 'SC with compact format should be empty';
};

for i = 1:size(test_cases, 1)
    grid_hole = test_cases{i, 1};
    brain_area = test_cases{i, 2};
    expected = test_cases{i, 3};
    description = test_cases{i, 4};

    result = parse_snc_subregion(grid_hole, brain_area);

    if strcmp(result, expected)
        fprintf('  [PASS] %s\n', description);
        test_passed = test_passed + 1;
    else
        fprintf('  [FAIL] %s\n', description);
        fprintf('         Expected: "%s", Got: "%s"\n', expected, result);
        test_failed = test_failed + 1;
    end
end

%% Test 2: Condition mask creation with actual session data
fprintf('\n--- Test 2: Condition Mask Creation ---\n');

% Load manifest
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
if ~exist(manifest_path, 'file')
    fprintf('  [SKIP] Manifest not found at: %s\n', manifest_path);
else
    manifest = readtable(manifest_path);

    % Find first SNc session with processed data
    snc_indices = find(strcmp(manifest.brain_area, 'SNc'));
    session_loaded = false;

    for idx = snc_indices'
        session_id = manifest.unique_id{idx};
        session_data_path = fullfile(project_root, 'data', 'processed', ...
            session_id, [session_id '_session_data.mat']);

        if exist(session_data_path, 'file')
            fprintf('  Loading session: %s\n', session_id);
            load(session_data_path, 'session_data');

            % Add metadata from manifest (as done in run_4factors_analysis.m)
            session_data.metadata = table2struct(manifest(idx, :));

            % Add SNc subregion to metadata
            session_data.metadata.snc_subregion = parse_snc_subregion(...
                session_data.metadata.grid_hole, ...
                session_data.metadata.brain_area);

            session_loaded = true;
            break;
        end
    end

    if ~session_loaded
        fprintf('  [SKIP] No processed SNc session data found\n');
    else
        % Call define_task_conditions
        try
            conditions = define_task_conditions(session_data);

            % Test: is_rvmSNc field exists
            if isfield(conditions, 'is_rvmSNc')
                fprintf('  [PASS] is_rvmSNc field exists\n');
                test_passed = test_passed + 1;
            else
                fprintf('  [FAIL] is_rvmSNc field missing\n');
                test_failed = test_failed + 1;
            end

            % Test: is_cdlSNc field exists
            if isfield(conditions, 'is_cdlSNc')
                fprintf('  [PASS] is_cdlSNc field exists\n');
                test_passed = test_passed + 1;
            else
                fprintf('  [FAIL] is_cdlSNc field missing\n');
                test_failed = test_failed + 1;
            end

            % Test: masks have correct length (same as other condition masks)
            n_trials = length(conditions.is_high_reward);
            if length(conditions.is_rvmSNc) == n_trials && ...
               length(conditions.is_cdlSNc) == n_trials
                fprintf('  [PASS] SNc masks have correct length (%d trials)\n', n_trials);
                test_passed = test_passed + 1;
            else
                fprintf('  [FAIL] SNc mask length mismatch\n');
                test_failed = test_failed + 1;
            end

            % Test: masks are mutually exclusive for SNc sessions
            expected_subregion = session_data.metadata.snc_subregion;
            if strcmp(expected_subregion, 'rvmSNc')
                if all(conditions.is_rvmSNc) && ~any(conditions.is_cdlSNc)
                    fprintf('  [PASS] rvmSNc session has correct mask values\n');
                    test_passed = test_passed + 1;
                else
                    fprintf('  [FAIL] rvmSNc session mask values incorrect\n');
                    test_failed = test_failed + 1;
                end
            elseif strcmp(expected_subregion, 'cdlSNc')
                if all(conditions.is_cdlSNc) && ~any(conditions.is_rvmSNc)
                    fprintf('  [PASS] cdlSNc session has correct mask values\n');
                    test_passed = test_passed + 1;
                else
                    fprintf('  [FAIL] cdlSNc session mask values incorrect\n');
                    test_failed = test_failed + 1;
                end
            end

        catch ME
            fprintf('  [FAIL] define_task_conditions threw error: %s\n', ME.message);
            test_failed = test_failed + 1;
        end
    end

    % Also test with an SC session if available
    sc_indices = find(strcmp(manifest.brain_area, 'SC'));
    sc_session_loaded = false;

    for idx = sc_indices'
        session_id = manifest.unique_id{idx};
        session_data_path = fullfile(project_root, 'data', 'processed', ...
            session_id, [session_id '_session_data.mat']);

        if exist(session_data_path, 'file')
            fprintf('\n  Loading SC session for comparison: %s\n', session_id);
            load(session_data_path, 'session_data');

            session_data.metadata = table2struct(manifest(idx, :));
            session_data.metadata.snc_subregion = parse_snc_subregion(...
                session_data.metadata.grid_hole, ...
                session_data.metadata.brain_area);

            sc_session_loaded = true;
            break;
        end
    end

    if sc_session_loaded
        try
            conditions = define_task_conditions(session_data);

            % Test: SC sessions should have both masks as false
            if ~any(conditions.is_rvmSNc) && ~any(conditions.is_cdlSNc)
                fprintf('  [PASS] SC session has both SNc masks as false\n');
                test_passed = test_passed + 1;
            else
                fprintf('  [FAIL] SC session should have both SNc masks as false\n');
                test_failed = test_failed + 1;
            end
        catch ME
            fprintf('  [FAIL] SC session test threw error: %s\n', ME.message);
            test_failed = test_failed + 1;
        end
    end
end

%% Test 3: Summary of all manifest sessions
fprintf('\n--- Test 3: SNc Subregion Assignment Summary ---\n');
fprintf('%-30s | %-6s | %-15s | %s\n', 'Session ID', 'Area', 'Grid Hole', 'Subregion');
fprintf('%s\n', repmat('-', 1, 75));

if exist('manifest', 'var')
    for i = 1:height(manifest)
        subregion = parse_snc_subregion(manifest.grid_hole{i}, manifest.brain_area{i});
        if isempty(subregion)
            subregion_display = '(n/a)';
        else
            subregion_display = subregion;
        end
        fprintf('%-30s | %-6s | %-15s | %s\n', ...
            manifest.unique_id{i}, manifest.brain_area{i}, ...
            manifest.grid_hole{i}, subregion_display);
    end
else
    fprintf('  Manifest not loaded - skipping summary\n');
end

%% Final Summary
fprintf('\n=== Test Summary ===\n');
fprintf('Passed: %d\n', test_passed);
fprintf('Failed: %d\n', test_failed);

if test_failed == 0
    fprintf('\nAll tests passed!\n');
else
    fprintf('\nSome tests failed. Please review the output above.\n');
end
