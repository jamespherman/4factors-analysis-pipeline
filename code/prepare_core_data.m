%% prepare_core_data.m
%
%   This function serves as the main wrapper for preparing analysis-ready
%   data structures for the 4factors task. It orchestrates the processing of
%   both neuronal and pupil data by calling specialized sub-functions.
%
% INPUTS:
%   session_data     - The main data structure containing all session info.
%   selected_neurons - A logical vector indicating which neurons to include.
%
% OUTPUT:
%   core_data        - A struct containing the processed neuronal and pupil
%                      data.
%
% Author: Jules
% Date: 2025-09-08

function core_data = prepare_core_data(session_data, selected_neurons, alignment_events)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Define Trial Indices
% The list of events to align to is now passed in as an argument.

% Load the structure containing all task codes
codes = initCodes();

% Find the indices of trials that belong to the '4factors' task. Importantly,
% we should only keep trials that are rewarded:
fourfactors_trial_indices = find(session_data.trialInfo.taskCode == ...
    codes.uniqueTaskCode_gSac_4factors & ~cellfun(@isempty, ...
    session_data.eventTimes.rewardCell) );

%% Prepare Neuronal Data
% Pass the alignment events to the neuronal data preparation function
core_data.spikes = prepare_neuronal_data(session_data, ...
    selected_neurons, fourfactors_trial_indices, alignment_events);

%% Prepare Pupil Data
% This part is removed as per the refactoring plan.
% core_data.pupil = prepare_pupil_data(session_data, ...
%     fourfactors_trial_indices, alignment_events);

end
