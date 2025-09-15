function [conditions, condition_defs] = define_task_conditions(varargin)
% DEFINE_TASK_CONDITIONS Creates trial condition masks and defines analysis plans
%
% This function serves a dual purpose for the 4factors task:
% 1.  When called without arguments, it returns a comprehensive analysis
%     plan in the `condition_defs` struct. This plan is the single source
%     of truth for all analyses run in the main pipeline.
% 2.  When called with session data, it calculates and returns a struct of
%     logical masks for various trial conditions based on the session's data.
%
% OUTPUTS:
%   conditions:     Struct of logical masks for trial conditions.
%   condition_defs: A struct containing the full, structured analysis plan.

%% --- I. Define the Comprehensive Analysis Plan ---
% This section defines the entire analysis plan. It is used by the main
% pipeline script, `run_4factors_analysis.m`, to determine which analyses to run.

% A. Alignment Events
condition_defs.alignment_events = {'fixOn', 'targetOn', 'fixOff', 'saccadeOnset', 'reward'};

% B. Canonical Names for Condition Masks
% These are the building blocks for the analysis plans below.
condition_defs.condition_masks.reward = {'is_high_reward', 'is_low_reward'};
condition_defs.condition_masks.salience = {'is_high_salience', 'is_low_salience'};
condition_defs.condition_masks.identity = {'is_face_target', 'is_nonface_target'};
condition_defs.condition_masks.probability = {'is_high_probability', 'is_low_probability'};
condition_defs.condition_masks.rf_location = {'is_in_rf', 'is_out_of_rf'};

% C. Baseline Comparison Plan
% Each element defines a "Baseline vs. Post-Event Activity" analysis.
baseline_conditions = [ ...
    condition_defs.condition_masks.reward, ...
    condition_defs.condition_masks.salience, ...
    condition_defs.condition_masks.identity, ...
    condition_defs.condition_masks.probability ...
];
condition_defs.baseline_plan = struct('name', {});
for i = 1:length(baseline_conditions)
    condition_defs.baseline_plan(i).name = baseline_conditions{i};
end

% D. ROC Analysis Plan
% Each element defines a bin-by-bin ROC comparison.
roc_plan_def = { ...
    'reward',      'targetOn', condition_defs.condition_masks.reward{1},      condition_defs.condition_masks.reward{2},      'is_in_rf'; ...
    'salience',    'targetOn', condition_defs.condition_masks.salience{1},    condition_defs.condition_masks.salience{2},    'is_in_rf'; ...
    'identity',    'targetOn', condition_defs.condition_masks.identity{1},    condition_defs.condition_masks.identity{2},    'is_in_rf'; ...
    'probability', 'targetOn', condition_defs.condition_masks.probability{1}, condition_defs.condition_masks.probability{2}, 'is_in_rf'  ...
};
condition_defs.roc_plan = struct('name', {}, 'event', {}, 'cond1', {}, 'cond2', {}, 'trial_mask', {});
for i = 1:size(roc_plan_def, 1)
    condition_defs.roc_plan(i).name       = roc_plan_def{i, 1};
    condition_defs.roc_plan(i).event      = roc_plan_def{i, 2};
    condition_defs.roc_plan(i).cond1      = roc_plan_def{i, 3};
    condition_defs.roc_plan(i).cond2      = roc_plan_def{i, 4};
    condition_defs.roc_plan(i).trial_mask = roc_plan_def{i, 5};
end

% E. Per-Neuron Diagnostic Plots
diag_plots(1).event = 'targetOn';
diag_plots(1).title = 'Target On';
diag_plots(1).conditions_to_compare = {'is_in_rf', 'is_out_of_rf'};

diag_plots(2).event = 'saccadeOnset';
diag_plots(2).title = 'Saccade Onset';
diag_plots(2).conditions_to_compare = {'is_in_rf', 'is_out_of_rf'};

condition_defs.diagnostic_plots = diag_plots;

% F. N-way ANOVA Plan
condition_defs.anova_plan.run = true;
condition_defs.anova_plan.event = 'targetOn';
condition_defs.anova_plan.factors = {'reward', 'salience', 'identity', 'probability'};
condition_defs.anova_plan.trial_mask = 'is_in_rf';

% G. Behavioral Analysis Plan
% Each element defines a behavioral analysis to be run.
condition_defs.behavior_plan(1).name = 'ReactionTime';
condition_defs.behavior_plan(1).dependent_variable = ...
    {'eventTimes.pdsSaccadeOn', '-', 'eventTimes.pdsFixOff'};
condition_defs.behavior_plan(1).factors = {'reward', 'salience', 'identity', 'probability'};
condition_defs.behavior_plan(1).trial_mask = 'is_in_rf';

condition_defs.behavior_plan(2).name = 'PeakVelocity';
condition_defs.behavior_plan(2).dependent_variable = {'trialInfo.peakVel'};
condition_defs.behavior_plan(2).factors = {'reward', 'salience', 'identity', 'probability'};
condition_defs.behavior_plan(2).trial_mask = 'is_in_rf';

condition_defs.behavior_plan(3).name = 'EndpointError';
% 'endpoint_error' is a special keyword that will be calculated by
% the analysis function from postSacXY and target location.
condition_defs.behavior_plan(3).dependent_variable = {'endpoint_error'};
condition_defs.behavior_plan(3).factors = {'reward', 'salience', 'identity', 'probability'};
condition_defs.behavior_plan(3).trial_mask = 'is_in_rf';


%% --- II. Mode Dispatch ---
% If the function is called without arguments, return the analysis plan.
if nargin == 0
    conditions = [];
    return;
end

%% --- III. Session-Specific Condition Mask Calculation ---
% This section runs only when session data is provided as input. It
% calculates the logical masks for each condition based on the trial data.
sessionData = varargin{1};
trialInfo = sessionData.trialInfo;
codes = initCodes;

% --- Initial Trial Filtering (The Master Mask) ---
% Identify trials belonging to the gSac_4factors task
isGSac4factors = trialInfo.taskCode == codes.uniqueTaskCode_gSac_4factors;


% Identify successfully completed trials (outcome == 1 is success)
isSuccessful = trialInfo.outcome == 1;

% The master mask: valid, completed trials from the correct task
masterMask = isGSac4factors & isSuccessful;

% --- Receptive Field (RF) Conditions ---
% Get the index of the target location inside the SC receptive field
rfLocationIdx = analysis.determine_rf_location(sessionData);

% Create logical masks based on the trial-by-trial target location
isInRf = (trialInfo.pdsTargLocIdx == rfLocationIdx);
isOutOfRf = ~isInRf;

% --- The Four Factors ---
% 1. Reward
isHighReward = (trialInfo.reward == 2);
isLowReward = (trialInfo.reward == 1);

% 2. Salience
isHighSalience = (trialInfo.salience == 2);
isLowSalience = (trialInfo.salience == 1);

% 3. Identity (stimType: 1=Face, 2=Non-face)
isFaceTarget = (trialInfo.stimType == 1);
isNonfaceTarget = (trialInfo.stimType == 2);

% 4. Spatial Probability (Data-Driven)
% Get target locations for all valid gSac_4factors trials
validTrialLocs = trialInfo.pdsTargLocIdx(masterMask);

% Calculate frequency for each unique target location
[uniqueLocs, ~, locIndices] = unique(validTrialLocs);
locCounts = accumarray(locIndices, 1);
medianFreq = median(locCounts);

% Identify high and low probability locations
highProbLocs = uniqueLocs(locCounts > medianFreq);
lowProbLocs = uniqueLocs(locCounts < medianFreq);

% Create masks based on whether trial's target is in a high/low prob
isHighProbability = ismember(trialInfo.pdsTargLocIdx, highProbLocs);
isLowProbability = ismember(trialInfo.pdsTargLocIdx, lowProbLocs);

% --- Output Structure ---
% Filter all logical masks by the masterMask and store in the output.
% This ensures all output fields have the same length, corresponding
% to the number of valid trials.
conditions.is_in_rf = isInRf(masterMask);
conditions.is_out_of_rf = isOutOfRf(masterMask);

conditions.is_high_reward = isHighReward(masterMask);
conditions.is_low_reward = isLowReward(masterMask);

conditions.is_high_salience = isHighSalience(masterMask);
conditions.is_low_salience = isLowSalience(masterMask);

conditions.is_face_target = isFaceTarget(masterMask);
conditions.is_nonface_target = isNonfaceTarget(masterMask);

conditions.is_high_probability = isHighProbability(masterMask);
conditions.is_low_probability = isLowProbability(masterMask);

end
