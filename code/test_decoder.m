%% test_decoder.m
%
% Tests a pre-trained SVM classifier on a specified test dataset.
%
% This function takes a collection of trained models and a testing plan
% item. It finds the correct model using the 'train_model_tag' from the
% testing plan, prepares the test data (features and labels), and then
% uses the model to predict labels for the test set. Finally, it
% calculates and returns the classification accuracy.
%
% ---
%
% Inputs:
%
%   trained_models: A cell array of modelInfo structs, where each struct
%                   is the output of train_decoder.m.
%
%   conditions: A struct with logical masks for experimental conditions.
%
%   core_data: A struct with processed neural data (e.g., firing rates).
%
%   testing_plan_item: A single struct from the analysis plan, defining
%                      the testing procedure. It must contain:
%                      - .train_model_tag: Identifier for the model to use.
%                      - .align_event: Event to align data to.
%                      - .time_window: [1x2] time window for features.
%                      - .cond1, .cond2: Conditions to decode.
%                      - .trial_mask: Mask for selecting trials.
%
% ---
%
% Output:
%
%   results: A struct containing the analysis results, including:
%            - .accuracy: The classification accuracy on the test set.
%
% ---
%
% Author: Jules
% Date: 2025-09-15

function results = test_decoder(trained_models, conditions, ...
    core_data, testing_plan_item)

%% Setup Paths
% Add the 'utils' directory to the path for helper functions.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Find Pre-Trained Model
% Find the correct model from the collection using the train_model_tag.
train_tag = testing_plan_item.train_model_tag;
model_found = false;
for i = 1:length(trained_models)
    if strcmp(trained_models{i}.model_tag, train_tag)
        modelInfo = trained_models{i};
        model_found = true;
        break;
    end
end

% Throw an error if the specified model was not found.
if ~model_found
    error('test_decoder:ModelNotFound', ...
        'Model with tag "%s" not found in trained_models.', train_tag);
end

%% Data Preparation
% Extract feature matrix (X_test) and label vector (Y_test).

% Select the binned firing rates for the correct alignment event
align_event = testing_plan_item.align_event;
binned_rates = core_data.(align_event).rates;
time_vector = core_data.(align_event).time_vector;

% Identify time bins within the specified time_window
time_window = testing_plan_item.time_window;
time_bin_indices = time_vector >= time_window(1) & ...
                   time_vector <= time_window(2);

% Average firing rates across time bins for the feature matrix
X = squeeze(mean(binned_rates(:, :, time_bin_indices), 3));

% Get logical masks for the two conditions and the trial mask
cond1_mask = conditions.(testing_plan_item.cond1);
cond2_mask = conditions.(testing_plan_item.cond2);
trial_mask = conditions.(testing_plan_item.trial_mask);

% Combine masks to select the final trials for testing
testing_mask = (cond1_mask | cond2_mask) & trial_mask;

% Select the feature matrix and create the label vector
X_test = X(testing_mask, :);
Y_test = cond1_mask(testing_mask);

%% Test Model
% Use the pre-trained model to predict labels and calculate accuracy.
Y_pred = predict(modelInfo.model, X_test);
accuracy = mean(Y_pred == Y_test);

%% Package Output
% Store the accuracy in the results struct.
results.accuracy = accuracy;

end
