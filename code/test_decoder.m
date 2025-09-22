%% test_decoder.m
%
% Tests a pre-trained SVM classifier based on a specified test plan.
%
% This function orchestrates different types of decoding tests based on the
% `testing_plan_item.type`. It supports cross-factor and cross-time
% generalization, as well as standard cross-validated accuracy estimation.
%
% ---
%
% Inputs:
%
%   trained_models: A cell array of modelInfo structs from train_decoder.m.
%
%   conditions: A struct with logical masks for experimental conditions.
%
%   core_data: A struct with processed neural data (e.g., firing rates).
%
%   testing_plan_item: A struct from the analysis plan defining the test.
%
% ---
%
% Output:
%
%   results: A struct with .accuracy and .accuracy_ci fields.
%
% ---
%
% Author: Jules
% Date: 2025-09-15

function results = test_decoder(trained_models, conditions, ...
    core_data, testing_plan_item)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Find Pre-Trained Model
train_tag = testing_plan_item.train_model_tag;
model_found = false;
for i = 1:length(trained_models)
    if strcmp(trained_models{i}.model_tag, train_tag)
        modelInfo = trained_models{i};
        model_found = true;
        break;
    end
end

if ~model_found
    error('test_decoder:ModelNotFound', ...
        'Model with tag "%s" not found.', train_tag);
end

%% Main Logic: Switch on Test Type
switch testing_plan_item.type
    case {'cross_factor', 'cross_time'}
        %% Data Preparation for Generalization Tests
        if strcmp(testing_plan_item.type, 'cross_time')
            event = testing_plan_item.test_event;
            time_window = testing_plan_item.test_time_window;
        else % 'cross_factor'
            event = testing_plan_item.event;
            time_window = testing_plan_item.time_window;
        end

        binned_rates = core_data.spikes.(event).rates;
        time_vector = core_data.spikes.(event).time_vector;

        time_bin_indices = time_vector >= time_window(1) & ...
                           time_vector <= time_window(2);
        X = squeeze(mean(binned_rates(:, :, time_bin_indices), 3));

        cond1_mask = conditions.(testing_plan_item.test_cond1);
        cond2_mask = conditions.(testing_plan_item.test_cond2);

        mask_field = testing_plan_item.trial_mask;
        if iscell(mask_field)
            trial_mask = true(size(cond1_mask));
            for i = 1:length(mask_field)
                trial_mask = trial_mask & conditions.(mask_field{i});
            end
        else
            trial_mask = conditions.(mask_field);
        end

        testing_mask = (cond1_mask | cond2_mask) & trial_mask;
        X_test = X(testing_mask, :);
        Y_test = cond1_mask(testing_mask);

        %% Test Model and Calculate Accuracy
        Y_pred = predict(modelInfo.model, X_test);
        n_correct = sum(Y_pred == Y_test);
        n_total = length(Y_test);
        [accuracy, accuracy_ci] = binofit(n_correct, n_total);

    case 'standard'
        %% Standard K-Fold Cross-Validation
        % Retrieve the original training data stored in the model object
        X_train = modelInfo.model.X;
        Y_train = modelInfo.model.Y;

        % Perform 10-fold cross-validation
        cv_model = crossval(modelInfo.model, 'KFold', 10);
        Y_pred = kfoldPredict(cv_model);

        % Calculate accuracy and confidence interval
        n_correct = sum(Y_pred == Y_train);
        n_total = length(Y_train);
        [accuracy, accuracy_ci] = binofit(n_correct, n_total);

    otherwise
        error('test_decoder:UnknownTestType', ...
            'Unknown test type: %s', testing_plan_item.type);
end

%% Package Output
results.accuracy = accuracy;
results.accuracy_ci = accuracy_ci;

end
