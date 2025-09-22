%% train_decoder.m
%
% Trains a suite of SVM classifiers based on a provided analysis plan.
%
% This function iterates through a set of training plan items, each
% defining a specific model to be trained. For each item, it prepares the
% feature matrix and label vector, then trains a ClassificationSVM model
% using the full dataset specified for that model. The trained model and
% its associated metadata are stored in a struct.
%
% ---
%
% Inputs:
%
%   session_data: The main data structure for a session.
%
%   conditions: A struct with logical masks for experimental conditions.
%
%   core_data: A struct with processed neural data (e.g., firing rates).
%
%   training_plan_item: A single struct from the analysis plan, defining
%                       the model to be trained. It must contain:
%                       - .model_tag: A unique identifier for the model.
%                       - .event: Event to align data to.
%                       - .time_window: [1x2] time window for features.
%                       - .cond1, .cond2: Conditions to decode.
%                       - .trial_mask: Mask for selecting trials.
%
% ---
%
% Output:
%
%   modelInfo: A struct containing the trained model and metadata:
%              - .model: The trained ClassificationSVM model object.
%              - .model_tag: The unique identifier for the model.
%              - (other metadata from training_plan_item is also stored)
%
% ---
%
% Author: Jules
% Date: 2025-09-15

function modelInfo = train_decoder(session_data, conditions, ...
    core_data, training_plan_item)

%% Setup Paths
% Add the 'utils' directory to the path for helper functions.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Data Preparation
% Extract feature matrix (X) and label vector (Y) for training.

% Select the binned firing rates for the correct alignment event
try
event = training_plan_item.event;
binned_rates = core_data.spikes.(event).rates;
time_vector = core_data.spikes.(event).time_vector;

% Identify time bins within the specified time_window
time_window = training_plan_item.time_window;
time_bin_indices = time_vector >= time_window(1) & ...
                   time_vector <= time_window(2);
catch me
    keyboard
end

% Average firing rates across time bins for the feature matrix
% The result is a [n_trials x n_neurons] matrix.
X = squeeze(mean(binned_rates(:, :, time_bin_indices), 3));

% Get logical masks for the two conditions
cond1_mask = conditions.(training_plan_item.cond1);
cond2_mask = conditions.(training_plan_item.cond2);

% Handle the trial_mask, which can be a single string or a cell array
mask_field = training_plan_item.trial_mask;
if iscell(mask_field)
    % If it's a cell, combine masks with logical AND
    trial_mask = true(size(cond1_mask)); % Start with all true
    for i = 1:length(mask_field)
        trial_mask = trial_mask & conditions.(mask_field{i});
    end
else
    % If it's a single string, retrieve the mask directly
    trial_mask = conditions.(mask_field);
end

% Combine masks to select the final trials for training
training_mask = (cond1_mask | cond2_mask) & trial_mask;

% Select the feature matrix and create the label vector
X_train = X(:, training_mask);
Y_train = cond1_mask(training_mask);

%% Train SVM Model
% Train a ClassificationSVM model on the full prepared dataset.
% A uniform prior is used to account for imbalanced trial counts.
svm_model = fitcsvm(X_train', Y_train(:), 'Prior', 'uniform');

%% Package Output
% Store the trained model and metadata in the output struct.
modelInfo = training_plan_item;
modelInfo.model = svm_model;

end
