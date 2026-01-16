%% test_window_roc.m
% Test script to verify window-based ROC analysis implementation
%
% This script loads a single session and runs analyze_window_roc to verify:
% 1. Output dimensions match neuron count
% 2. AUC values bounded [0, 1]
% 3. P-values bounded [0, 1]
% 4. Salience results are NaN for non-bullseye sessions (if applicable)
% 5. Identity results are NaN for non-image sessions (if applicable)
% 6. Confidence intervals bracket AUC values

clear; clc; close all;

%% Setup paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

%% Select test session
% Test with one SC session and one SNc session
test_sessions = {'Feynman_08_12_2025_SC', 'Feynman_08_12_2025_SNc'};

%% Load analysis plan
fprintf('Loading analysis plan...\n');
[~, analysis_plan] = define_task_conditions();
fprintf('Analysis plan loaded.\n\n');

%% Run tests for each session
for i_session = 1:length(test_sessions)
    session_id = test_sessions{i_session};
    fprintf('==============================================\n');
    fprintf('Testing session: %s\n', session_id);
    fprintf('==============================================\n\n');

    %% Load session data
    local_processed_dir = fullfile(project_root, 'data', 'processed', session_id);
    local_session_data_path = fullfile(local_processed_dir, [session_id, '_session_data.mat']);

    if ~exist(local_session_data_path, 'file')
        fprintf('Session data not found at: %s\n', local_session_data_path);
        fprintf('Trying OneDrive location...\n');
        one_drive_path = findOneDrive;
        source_path = fullfile(one_drive_path, 'Neuronal Data Analysis', ...
            session_id, [session_id '_session_data.mat']);
        if ~exist(source_path, 'file')
            fprintf('SKIP: Session data not found.\n\n');
            continue;
        end
        load(source_path, 'session_data');
    else
        load(local_session_data_path, 'session_data');
    end
    fprintf('Session data loaded.\n');

    %% Get metadata from manifest
    manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
    manifest = readtable(manifest_path);
    session_idx = find(strcmp(manifest.unique_id, session_id));
    session_data.metadata = table2struct(manifest(session_idx, :));
    session_data.metadata.snc_subregion = parse_snc_subregion(...
        session_data.metadata.grid_hole, session_data.metadata.brain_area);

    fprintf('Brain area: %s\n', session_data.metadata.brain_area);

    %% Define conditions
    fprintf('Defining task conditions...\n');
    conditions = define_task_conditions(session_data);
    fprintf('Conditions defined.\n');

    %% Check if core_data exists, otherwise skip
    if ~isfield(session_data, 'analysis') || ~isfield(session_data.analysis, 'core_data')
        fprintf('SKIP: core_data not found. Run main pipeline first.\n\n');
        continue;
    end
    core_data = session_data.analysis.core_data;

    %% Get number of neurons
    first_event = analysis_plan.window_roc_plan.epoch_windows.visual.event;
    n_neurons = size(core_data.spikes.(first_event).rates, 1);
    fprintf('Number of neurons: %d\n\n', n_neurons);

    %% Run window ROC analysis
    fprintf('Running analyze_window_roc...\n');
    tic;
    window_roc = analyze_window_roc(session_data, conditions, analysis_plan, core_data);
    elapsed = toc;
    fprintf('Analysis completed in %.2f seconds.\n\n', elapsed);

    %% Verify outputs
    fprintf('--- VERIFICATION RESULTS ---\n\n');

    epoch_names = fieldnames(window_roc);
    all_tests_passed = true;

    for i_epoch = 1:length(epoch_names)
        epoch_name = epoch_names{i_epoch};
        fprintf('Epoch: %s\n', epoch_name);

        factor_names = fieldnames(window_roc.(epoch_name));
        for i_factor = 1:length(factor_names)
            factor_name = factor_names{i_factor};
            result = window_roc.(epoch_name).(factor_name);

            fprintf('  Factor: %s\n', factor_name);

            % 1. Check dimensions
            dim_check = size(result.auc, 1) == n_neurons && ...
                        size(result.p, 1) == n_neurons && ...
                        size(result.ci, 1) == n_neurons && ...
                        size(result.n_trials, 1) == n_neurons;
            if dim_check
                fprintf('    [PASS] Dimensions match neuron count (%d)\n', n_neurons);
            else
                fprintf('    [FAIL] Dimension mismatch!\n');
                all_tests_passed = false;
            end

            % 2. Check AUC bounds [0, 1]
            valid_auc = result.auc(~isnan(result.auc));
            if isempty(valid_auc)
                fprintf('    [WARN] All AUC values are NaN\n');
            else
                auc_in_bounds = all(valid_auc >= 0 & valid_auc <= 1);
                if auc_in_bounds
                    fprintf('    [PASS] AUC values in [0,1] (range: [%.3f, %.3f])\n', ...
                        min(valid_auc), max(valid_auc));
                else
                    fprintf('    [FAIL] AUC values out of bounds!\n');
                    all_tests_passed = false;
                end
            end

            % 3. Check p-value bounds [0, 1]
            valid_p = result.p(~isnan(result.p));
            if isempty(valid_p)
                fprintf('    [WARN] All p-values are NaN\n');
            else
                p_in_bounds = all(valid_p >= 0 & valid_p <= 1);
                if p_in_bounds
                    fprintf('    [PASS] P-values in [0,1] (range: [%.3f, %.3f])\n', ...
                        min(valid_p), max(valid_p));
                else
                    fprintf('    [FAIL] P-values out of bounds!\n');
                    all_tests_passed = false;
                end
            end

            % 4. Check CI brackets AUC
            valid_idx = ~isnan(result.auc) & ~isnan(result.ci(:,1)) & ~isnan(result.ci(:,2));
            if any(valid_idx)
                ci_brackets = all(result.ci(valid_idx,1) <= result.auc(valid_idx) & ...
                                  result.auc(valid_idx) <= result.ci(valid_idx,2));
                if ci_brackets
                    fprintf('    [PASS] CI brackets AUC values\n');
                else
                    fprintf('    [FAIL] CI does not bracket AUC!\n');
                    all_tests_passed = false;
                end
            else
                fprintf('    [WARN] No valid CI values to check\n');
            end

            % 5. Report trial counts and valid neuron count
            n_valid = sum(~isnan(result.auc));
            mean_trials = mean(result.n_trials(~isnan(result.auc), :), 1, 'omitnan');
            fprintf('    [INFO] Valid neurons: %d/%d, Mean trials: [%.1f, %.1f]\n', ...
                n_valid, n_neurons, mean_trials(1), mean_trials(2));

            % 6. Check factor-specific subset behavior
            if strcmp(factor_name, 'salience')
                n_bullseye = sum(conditions.is_bullseye_target);
                fprintf('    [INFO] Bullseye trials: %d\n', n_bullseye);
                if n_bullseye == 0 && all(isnan(result.auc))
                    fprintf('    [PASS] Correctly returns NaN for non-bullseye session\n');
                end
            elseif strcmp(factor_name, 'identity')
                n_image = sum(conditions.is_image_target);
                fprintf('    [INFO] Image trials: %d\n', n_image);
                if n_image == 0 && all(isnan(result.auc))
                    fprintf('    [PASS] Correctly returns NaN for non-image session\n');
                end
            end
        end
        fprintf('\n');
    end

    %% Summary
    fprintf('--- SUMMARY ---\n');
    if all_tests_passed
        fprintf('[PASS] All verification tests passed for %s!\n\n', session_id);
    else
        fprintf('[FAIL] Some tests failed for %s. Review output above.\n\n', session_id);
    end

    %% Display sample output structure
    fprintf('--- SAMPLE OUTPUT STRUCTURE ---\n');
    fprintf('window_roc fields: %s\n', strjoin(epoch_names, ', '));
    fprintf('window_roc.visual fields: %s\n', strjoin(fieldnames(window_roc.visual), ', '));
    fprintf('window_roc.visual.reward fields: %s\n', ...
        strjoin(fieldnames(window_roc.visual.reward), ', '));
    fprintf('\n');
end

fprintf('==============================================\n');
fprintf('Test script completed.\n');
fprintf('==============================================\n');
