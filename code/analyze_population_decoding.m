%% analyze_population_decoding.m
%
% This function performs population decoding analysis using a support vector
% machine (SVM) classifier. It is designed to be called by a main analysis
% script and operates on a single decoding plan item. The function supports
% three types of decoding: standard, cross-factor, and cross-time.
%
% The feature for each neuron on a given trial is the mean firing rate
% across the time bins defined by the epoch's time window. The SVM is
% always trained with a uniform prior to account for imbalanced trial
% counts between conditions.
%
% ---
%
% Inputs:
%
%   session_data: The main data structure for a session, typically loaded
%                 from a session_data.mat file.
%
%   conditions: A struct containing logical masks for different experimental
%               conditions (e.g., `conditions.is_high_reward`). These masks
%               are used to select trials for training and testing.
%
%   core_data: A struct containing the processed neural data, including
%              binned firing rates (`core_data.eventName.rates`) and their
%              corresponding time vectors (`core_data.eventName.time_vector`).
%
%   plan_item: A struct that specifies the parameters for the decoding
%              analysis. It must contain the following fields:
%              - .type: 'standard', 'cross_factor', or 'cross_time'
%              - .align_event: The event to align the analysis to.
%              - .time_window: A [1x2] vector defining the start and end
%                              of the time window for feature calculation.
%              - .cond1, .cond2: The names of the two conditions to be
%                                decoded (e.g., 'is_high_reward').
%              - .trial_mask: A mask to select a subset of trials for the
%                             analysis.
%              - .train_factor, .test_factor (for 'cross_factor' type)
%              - .train_epoch, .test_epoch (for 'cross_time' type)
%
% ---
%
% Output:
%
%   results: A struct containing the results of the decoding analysis.
%            The primary field is:
%            - .accuracy: The classification accuracy, which can be either
%                         cross-validated or based on a separate test set,
%                         depending on the analysis type.
%
% ---
%
% Author: Jules
% Date: 2025-09-15

function results = analyze_population_decoding(session_data, conditions, ...
    core_data, plan_item)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Data Preparation
% This section handles the initial data preparation which is common to
% most analysis types. The 'cross_time' analysis is a special case and
% will call a helper function to prepare its data since it needs to
% perform this process twice.

% A helper function to extract features for a given epoch. This is
% defined here to be accessible to all parts of the main function.
function X_epoch = extract_features(core_data, align_event, time_window)
    % Select the binned firing rates for the correct alignment event
    binned_rates = core_data.(align_event).rates;
    time_vector = core_data.(align_event).time_vector;

    % Identify the time bins that fall within the specified time_window
    time_bin_indices = time_vector >= time_window(1) & ...
                       time_vector <= time_window(2);

    % Average the firing rates across these time bins to create a single
    % feature per neuron for each trial.
    % The result is a [n_trials x n_neurons] matrix.
    X_epoch = squeeze(mean(binned_rates(:, :, time_bin_indices), 3));
end

% For standard and cross-factor, we prepare one primary feature matrix.
% The cross_time case will handle its own feature extraction.
if ~strcmp(plan_item.type, 'cross_time')
    X = extract_features(core_data, plan_item.align_event, ...
        plan_item.time_window);
end


%% Main Analysis Logic
switch plan_item.type
    case 'standard'
        % --- Standard Decoding: Train and test on the same dataset ---

        % Get the logical masks for the two conditions and the trial mask
        cond1_mask = conditions.(plan_item.cond1);
        cond2_mask = conditions.(plan_item.cond2);
        trial_mask = conditions.(plan_item.trial_mask);

        % Combine masks to select the final trials for classification
        final_trial_mask = (cond1_mask | cond2_mask) & trial_mask;

        % Select the feature matrix and create the label vector
        X_selected = X(final_trial_mask, :);
        Y = cond1_mask(final_trial_mask);

        % Train the SVM classifier with 10-fold cross-validation
        % A uniform prior is used to account for imbalanced trial counts.
        svm_model = fitcsvm(X_selected, Y, 'CrossVal', 'on', ...
            'KFold', 10, 'Prior', 'uniform');

        % Calculate the cross-validated classification accuracy
        results.accuracy = 1 - kfoldLoss(svm_model);
    case 'cross_factor'
        % --- Cross-Factor Decoding: Train on one factor, test on another ---

        % Get masks for conditions and overall trial mask
        cond1_mask = conditions.(plan_item.cond1);
        cond2_mask = conditions.(plan_item.cond2);
        trial_mask = conditions.(plan_item.trial_mask);

        % Get masks for training and testing factors
        train_factor_mask = conditions.(plan_item.train_factor);
        test_factor_mask = conditions.(plan_item.test_factor);

        % Define training trials
        train_mask = (cond1_mask | cond2_mask) & trial_mask & train_factor_mask;
        X_train = X(train_mask, :);
        Y_train = cond1_mask(train_mask);

        % Define testing trials
        test_mask = (cond1_mask | cond2_mask) & trial_mask & test_factor_mask;
        X_test = X(test_mask, :);
        Y_test = cond1_mask(test_mask);

        % Train the SVM classifier on the training set
        svm_model = fitcsvm(X_train, Y_train, 'Prior', 'uniform');

        % Predict labels for the testing set
        Y_pred = predict(svm_model, X_test);

        % Calculate accuracy
        results.accuracy = mean(Y_pred == Y_test);
    case 'cross_time'
        % --- Cross-Time Decoding: Train on one epoch, test on another ---

        % Extract features for the training epoch
        X_train_epoch = extract_features(core_data, ...
            plan_item.train_epoch.align_event, ...
            plan_item.train_epoch.time_window);

        % Extract features for the testing epoch
        X_test_epoch = extract_features(core_data, ...
            plan_item.test_epoch.align_event, ...
            plan_item.test_epoch.time_window);

        % Get masks for conditions and overall trial mask
        cond1_mask = conditions.(plan_item.cond1);
        cond2_mask = conditions.(plan_item.cond2);
        trial_mask = conditions.(plan_item.trial_mask);

        % Define the trials to be used for classification
        final_trial_mask = (cond1_mask | cond2_mask) & trial_mask;

        % Select training and testing data
        X_train = X_train_epoch(final_trial_mask, :);
        X_test = X_test_epoch(final_trial_mask, :);
        Y_labels = cond1_mask(final_trial_mask);

        % Train the SVM classifier on the training epoch data
        svm_model = fitcsvm(X_train, Y_labels, 'Prior', 'uniform');

        % Predict labels for the testing epoch data
        Y_pred = predict(svm_model, X_test);

        % Calculate accuracy
        results.accuracy = mean(Y_pred == Y_labels);
end

end
