function conditions = define_task_conditions(session_data)
% DEFINE_TASK_CONDITIONS creates a struct of logical masks for experimental conditions.
%
% This function identifies all valid, successfully completed trials from the
% gSac_4factors task and generates logical masks for various experimental
% conditions. All output masks are filtered by this master set of valid trials.
%
% INPUT:
%   session_data: A struct containing session data, including trialInfo.
%
% OUTPUT:
%   conditions: A struct with fields for each logical mask. Each field will
%               have a length equal to the number of valid gSac_4factors trials.
%
% See also: analysis.determine_rf_location

% --- Initial Trial Filtering (The Master Mask) ---
trialInfo = session_data.trialInfo;

% Identify trials belonging to the gSac_4factors task
is_gSac_4factors = strcmp(trialInfo.task, 'gSac_4factors');

% Identify successfully completed trials (assuming outcome == 1 is success)
is_successful = trialInfo.outcome == 1;

% The master mask: valid, completed trials from the correct task
master_mask = is_gSac_4factors & is_successful;

% --- Receptive Field (RF) Conditions ---
% Get the index of the target location inside the SC receptive field
rf_location_idx = analysis.determine_rf_location(session_data);

% Create logical masks based on the trial-by-trial target location
is_in_rf = (trialInfo.pdsTargLocIdx == rf_location_idx);
is_out_of_rf = ~is_in_rf;

% --- The Four Factors ---

% 1. Reward
is_high_reward = (trialInfo.reward == 2);
is_low_reward = (trialInfo.reward == 1);

% 2. Salience
is_high_salience = (trialInfo.salience == 2);
is_low_salience = (trialInfo.salience == 1);

% 3. Identity (stimType: 1=Face, 2=Non-face)
is_face_target = (trialInfo.stimType == 1);
is_nonface_target = (trialInfo.stimType == 2);

% 4. Spatial Probability (Data-Driven)
% Get target locations for all valid gSac_4factors trials
valid_trial_locs = trialInfo.pdsTargLocIdx(master_mask);

% Calculate frequency of presentation for each unique target location
[unique_locs, ~, loc_indices] = unique(valid_trial_locs);
loc_counts = accumarray(loc_indices, 1);
median_freq = median(loc_counts);

% Identify high and low probability locations
high_prob_locs = unique_locs(loc_counts > median_freq);
low_prob_locs = unique_locs(loc_counts < median_freq);

% Create masks based on whether the trial's target is in a high/low prob location
is_high_probability = ismember(trialInfo.pdsTargLocIdx, high_prob_locs);
is_low_probability = ismember(trialInfo.pdsTargLocIdx, low_prob_locs);

% --- Output Structure ---
% Filter all logical masks by the master_mask and store in the output struct.
% This ensures all output fields have the same length, corresponding to the
% number of valid trials.
conditions.is_in_rf = is_in_rf(master_mask);
conditions.is_out_of_rf = is_out_of_rf(master_mask);

conditions.is_high_reward = is_high_reward(master_mask);
conditions.is_low_reward = is_low_reward(master_mask);

conditions.is_high_salience = is_high_salience(master_mask);
conditions.is_low_salience = is_low_salience(master_mask);

conditions.is_face_target = is_face_target(master_mask);
conditions.is_nonface_target = is_nonface_target(master_mask);

conditions.is_high_probability = is_high_probability(master_mask);
conditions.is_low_probability = is_low_probability(master_mask);

end
