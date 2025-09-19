function [conditions, condition_defs] = define_task_conditions(varargin)
% DEFINE_TASK_CONDITIONS Creates trial condition masks and defines
% analysis plans
%
% This function serves a dual purpose for the 4factors task:
% 1.  When called without arguments, it returns a comprehensive
%     analysis plan in the `condition_defs` struct. This plan is
%     the single source of truth for all analyses run in the main
%     pipeline.
% 2.  When called with session data, it determines the recorded
%     hemisphere based on manifest data and calculates a struct of
%     logical masks for various spatial and trial-based conditions.
%
% OUTPUTS:
%   conditions:     A struct of logical masks for trial conditions,
%                   specific to a single session's data. Each field
%                   is a logical vector where `true` indicates that
%                   a trial meets a specific condition (e.g.,
%                   `is_high_reward`). These masks are filtered to
%                   only include valid, successful trials from the
%                   '4factors' task.
%
%   condition_defs: A struct containing the full, session-agnostic
%                   analysis_plan. This serves as the single source
%                   of truth for what analyses to run, what
%                   conditions to compare, and how to configure
%                   them. It is used by the main analysis pipeline
%                   to orchestrate the entire analysis. The
%                   `.decoding_plan` field is further divided into a
%                   `.training_plan` for defining classifiers and a
%                   `.testing_plan` for defining generalization
%                   tests.

%% --- I. Define the Comprehensive Analysis Plan ---
% This section defines the entire analysis plan. It is used by the
% main pipeline script, `run_4factors_analysis.m`, to determine
% which analyses to run.

% A. Alignment Events
condition_defs.alignment_events = {'fixOn', 'targetOn', 'fixOff', ...
    'saccadeOnset', 'reward'};

% B. Canonical Names for Condition Masks
% These are the building blocks for the analysis plans below.
condition_defs.condition_masks.reward = ...
    {'is_high_reward', 'is_low_reward'};
condition_defs.condition_masks.salience = ...
    {'is_high_salience', 'is_low_salience'};
condition_defs.condition_masks.identity = ...
    {'is_face_target', 'is_nonface_target'};
condition_defs.condition_masks.probability = ...
    {'is_high_probability', 'is_low_probability'};
condition_defs.condition_masks.spatial = {'is_contralateral_target', ...
    'is_ipsilateral_target', 'is_opposite_rf'};

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
    'reward', 'targetOn', ...
    condition_defs.condition_masks.reward{1}, ...
    condition_defs.condition_masks.reward{2}, ...
    'is_contralateral_target'; ...
    'salience', 'targetOn', ...
    condition_defs.condition_masks.salience{1}, ...
    condition_defs.condition_masks.salience{2}, ...
    'is_contralateral_target'; ...
    'identity', 'targetOn', ...
    condition_defs.condition_masks.identity{1}, ...
    condition_defs.condition_masks.identity{2}, ...
    'is_contralateral_target'; ...
    'probability', 'targetOn', ...
    condition_defs.condition_masks.probability{1}, ...
    condition_defs.condition_masks.probability{2}, ...
    'is_contralateral_target' ...
};
condition_defs.roc_plan = struct('name', {}, 'event', {}, ...
    'cond1', {}, 'cond2', {}, 'trial_mask', {});
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
diag_plots(1).conditions_to_compare = ...
    {'is_contralateral_target', 'is_ipsilateral_target'};
diag_plots(2).event = 'saccadeOnset';
diag_plots(2).title = 'Saccade Onset';
diag_plots(2).conditions_to_compare = ...
    {'is_contralateral_target', 'is_ipsilateral_target'};
condition_defs.diagnostic_plots = diag_plots;

% F. N-way ANOVA Plan
condition_defs.anova_plan.run = true;
condition_defs.anova_plan.event = 'targetOn';
condition_defs.anova_plan.factors = ...
    {'reward', 'salience', 'identity', 'probability'};
condition_defs.anova_plan.trial_mask = 'is_contralateral_target';

% G. Behavioral Analysis Plan
% Each element defines a behavioral analysis to be run.
condition_defs.behavior_plan(1).name = 'ReactionTime';
condition_defs.behavior_plan(1).dependent_variable = ...
    {'eventTimes.pdsSaccadeOn', '-', 'eventTimes.pdsFixOff'};
condition_defs.behavior_plan(1).factors = ...
    {'reward', 'salience', 'identity', 'probability'};
condition_defs.behavior_plan(1).trial_mask = 'is_contralateral_target';
condition_defs.behavior_plan(2).name = 'PeakVelocity';
condition_defs.behavior_plan(2).dependent_variable = ...
    {'trialInfo.peakVel'};
condition_defs.behavior_plan(2).factors = ...
    {'reward', 'salience', 'identity', 'probability'};
condition_defs.behavior_plan(2).trial_mask = 'is_contralateral_target';
condition_defs.behavior_plan(3).name = 'EndpointError';
% 'endpoint_error' is a special keyword that will be calculated by
% the analysis function from postSacXY and target location.
condition_defs.behavior_plan(3).dependent_variable = {'endpoint_error'};
condition_defs.behavior_plan(3).factors = ...
    {'reward', 'salience', 'identity', 'probability'};
condition_defs.behavior_plan(3).trial_mask = 'is_contralateral_target';

% H. Population Decoding Plan
% This plan is now a two-stage process. First, a set of classifiers
% are trained based on the 'training_plan'. Second, these trained
% models are tested on new data according to the 'testing_plan'.

% H1. Training Plan
% Defines all classifiers to be trained.
training_plan = struct('model_tag', {}, 'align_event', {}, ...
                       'time_window', {}, 'cond1', {}, 'cond2', {}, ...
                       'trial_mask', {});
factors = {'reward', 'salience', 'identity', 'probability', 'spatial'};
epoch_defs = struct(...
    'Visual', struct('event', 'targetOn', 'window', [0.05, 0.25]), ...
    'Delay',  struct('event', 'fixOff', 'window', [-0.15, 0.05]), ...
    'Saccade',struct('event', 'saccadeOnset', 'window', [-0.20, 0.00]) ...
);
epoch_names = fieldnames(epoch_defs);
for i = 1:length(factors)
    factor_name = factors{i};
    conds = condition_defs.condition_masks.(factor_name);
    for j = 1:length(epoch_names)
        epoch_name = epoch_names{j};
        epoch_details = epoch_defs.(epoch_name);
        entry = struct();
        % Capitalize first letter for tag
        tag_factor_name = [upper(factor_name(1)), factor_name(2:end)];
        if strcmp(factor_name, 'spatial')
            tag_factor_name = 'Spatial';
        end
        entry.model_tag = sprintf('%s_%s', tag_factor_name, epoch_name);
        entry.align_event = epoch_details.event;
        entry.time_window = epoch_details.window;
        entry.cond1 = conds{1};
        entry.cond2 = conds{2};
        if strcmp(factor_name, 'identity')
            entry.trial_mask = ...
                {'is_contralateral_target', 'is_image_target'};
        else
            entry.trial_mask = 'is_contralateral_target';
        end
        training_plan(end+1) = entry;
    end
end
condition_defs.decoding_plan.training_plan = training_plan;

% H2. Testing Plan
% Defines all generalization tests to be performed.
testing_plan = struct('test_name', {}, 'type', {}, ...
                      'train_model_tag', {}, 'test_cond1', {}, ...
                      'test_cond2', {}, 'align_event', {}, ...
                      'time_window', {}, 'test_align_event', {}, ...
                      'test_time_window', {}, 'trial_mask', {});

% --- Cross-Factor Tests ---
cross_factors = {'Reward', 'Salience', 'Identity', 'Probability'};
cross_epochs = {'Visual', 'Saccade'};
for i = 1:length(cross_factors)
    train_factor = cross_factors{i};
    for j = 1:length(cross_factors)
        if i == j, continue; end
        test_factor = cross_factors{j};
        for k = 1:length(cross_epochs)
            epoch_name = cross_epochs{k};
            epoch_details = epoch_defs.(epoch_name);
            entry = struct();
            entry.test_name = sprintf('%s_x_%s_%s', train_factor, ...
                                      test_factor, epoch_name);
            entry.type = 'cross_factor';
            entry.train_model_tag = sprintf('%s_%s', train_factor, ...
                                            epoch_name);
            test_conds = ...
                condition_defs.condition_masks.(lower(test_factor));
            entry.test_cond1 = test_conds{1};
            entry.test_cond2 = test_conds{2};
            entry.align_event = epoch_details.event;
            entry.time_window = epoch_details.window;
            % Fields for cross-time, leave empty
            entry.test_align_event = '';
            entry.test_time_window = [];
            if strcmp(train_factor, 'Identity') || ...
               strcmp(test_factor, 'Identity')
                entry.trial_mask = ...
                    {'is_contralateral_target', 'is_image_target'};
            else
                entry.trial_mask = 'is_contralateral_target';
            end
            testing_plan(end+1) = entry;
        end
    end
end

% --- Cross-Time Tests ---
main_factors = {'Reward', 'Salience', 'Identity', 'Probability'};
for i = 1:length(main_factors)
    factor = main_factors{i};
    factor_conds = condition_defs.condition_masks.(lower(factor));
    % Visual -> Saccade
    entry = struct();
    entry.test_name = sprintf('%s_Vis_x_Sac', factor);
    entry.type = 'cross_time';
    entry.train_model_tag = sprintf('%s_Visual', factor);
    entry.test_align_event = epoch_defs.Saccade.event;
    entry.test_time_window = epoch_defs.Saccade.window;
    entry.test_cond1 = factor_conds{1};
    entry.test_cond2 = factor_conds{2};
    % Fields for cross-factor, leave empty
    entry.align_event = '';
    entry.time_window = [];
    if strcmp(factor, 'Identity')
        entry.trial_mask = ...
            {'is_contralateral_target', 'is_image_target'};
    else
        entry.trial_mask = 'is_contralateral_target';
    end
    testing_plan(end+1) = entry;
    % Saccade -> Visual
    entry = struct();
    entry.test_name = sprintf('%s_Sac_x_Vis', factor);
    entry.type = 'cross_time';
    entry.train_model_tag = sprintf('%s_Saccade', factor);
    entry.test_align_event = epoch_defs.Visual.event;
    entry.test_time_window = epoch_defs.Visual.window;
    entry.test_cond1 = factor_conds{1};
    entry.test_cond2 = factor_conds{2};
    % Fields for cross-factor, leave empty
    entry.align_event = '';
    entry.time_window = [];
    if strcmp(factor, 'Identity')
        entry.trial_mask = ...
            {'is_contralateral_target', 'is_image_target'};
    else
        entry.trial_mask = 'is_contralateral_target';
    end
    testing_plan(end+1) = entry;
end
condition_defs.decoding_plan.testing_plan = testing_plan;

% --- Standard CV Tests ---
% For every trained model, create a 'standard' test for it.
for i = 1:length(training_plan)
    train_item = training_plan(i);
    entry = struct();
    entry.test_name = sprintf('%s_CV', train_item.model_tag);
    entry.type = 'standard';
    entry.train_model_tag = train_item.model_tag;
    % --- Fill in dummy values for unused fields ---
    entry.test_cond1 = '';
    entry.test_cond2 = '';
    entry.align_event = '';
    entry.time_window = [];
    entry.test_align_event = '';
    entry.test_time_window = [];
    entry.trial_mask = '';
    testing_plan(end+1) = entry;
end
condition_defs.decoding_plan.testing_plan = testing_plan;

%% --- II. Mode Dispatch ---
% If the function is called without arguments, return the analysis
% plan.
if nargin == 0
    conditions = [];
    return;
end

%% --- III. Session-Specific Condition Mask Calculation ---
% This section runs only when session data is provided as input. It
% calculates the logical masks for each condition based on the trial
% data.
sessionData = varargin{1};
trialInfo = sessionData.trialInfo;
eventTimes = sessionData.eventTimes;
codes = initCodes;

% --- Initial Trial Filtering (The Master Mask) ---
% Identify trials belonging to the gSac_4factors task
isGSac4factors = ...
    trialInfo.taskCode == codes.uniqueTaskCode_gSac_4factors;
% Identify successfully completed trials (outcome == 1 is success)
isSuccessful = ~cellfun(@isempty, eventTimes.rewardCell);
% The master mask: valid, completed trials from the correct task
masterMask = isGSac4factors & isSuccessful;

% --- Spatial Conditions from Ground Truth ---
% Determine recorded hemisphere by parsing the 'grid_hole' string
% from the session manifest. This provides the ground truth for the
% recording location.
grid_hole_str = sessionData.metadata.grid_hole;

% Robustly parse the (X, Y) coordinates from the string. The sscanf
% function is used here to extract the first floating-point number,
% which corresponds to the X-coordinate. It is designed to handle
% variations in spacing and the presence of parentheses and commas.
coords = sscanf(grid_hole_str, '(%f, %f)');
if ~isempty(coords)
    grid_x = coords(1);
else
    grid_x = 0; % Default if parsing fails
    warning('Could not parse grid_hole string: "%s". Defaulting grid_x to 0.', ...
        grid_hole_str);
end

if grid_x < 0
    scSide = 'left';
    % Left SC -> Contralateral is Right Visual Field (theta > 0)
    contra_thetas = trialInfo.targetTheta > 0;
elseif grid_x > 0
    scSide = 'right';
    % Right SC -> Contralateral is Left Visual Field (theta < 0)
    contra_thetas = trialInfo.targetTheta < 0;
else
    % Default to right visual field if grid_x is 0 or ambiguous
    scSide = 'unknown';
    contra_thetas = trialInfo.targetTheta > 0;
    warning(['grid_x is 0 or could not be determined; defaulting ' ...
        'contralateral to right visual field.']);
end
isContralateral = contra_thetas;
isIpsilateral = ~contra_thetas;

% --- Opposite RF Mask ---
% Find the most frequent target location in the contralateral
% hemifield
validContraTrials = masterMask & isContralateral;
contraLocations = trialInfo.targetLocIdx(validContraTrials);
if ~isempty(contraLocations)
    primaryContraLocation = mode(contraLocations);
    % Determine the location opposite to the primary one
    switch primaryContraLocation
        case 1
            oppositeLocation = 3;
        case 3
            oppositeLocation = 1;
        case 2
            oppositeLocation = 4;
        case 4
            oppositeLocation = 2;
        otherwise
            oppositeLocation = NaN; % Should not happen
    end
    % Create the mask for trials at the opposite location
    isOppositeRf = (trialInfo.targetLocIdx == oppositeLocation);
else
    % If no contralateral trials, the mask is all false
    isOppositeRf = false(size(trialInfo.targetLocIdx));
end

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
is_image_target = isFaceTarget | isNonfaceTarget;

% 4. Spatial Probability (Data-Driven)
% Get target locations for all valid gSac_4factors trials
validTrialLocs = trialInfo.targetLocIdx(masterMask);
% Calculate frequency for each unique target location
[uniqueLocs, ~, locIndices] = unique(validTrialLocs);
locCounts = accumarray(locIndices, 1);
medianFreq = median(locCounts);
% Identify high and low probability locations
highProbLocs = uniqueLocs(locCounts > medianFreq);
lowProbLocs = uniqueLocs(locCounts < medianFreq);
% Create masks based on whether trial's target is in a high/low prob
isHighProbability = ismember(trialInfo.targetLocIdx, highProbLocs);
isLowProbability = ismember(trialInfo.targetLocIdx, lowProbLocs);

% --- Output Structure ---
% Filter all logical masks by the masterMask and store in the output.
% This ensures all output fields have the same length, corresponding
% to the number of valid trials.
conditions.is_contralateral_target = isContralateral(masterMask);
conditions.is_ipsilateral_target = isIpsilateral(masterMask);
conditions.is_opposite_rf = isOppositeRf(masterMask);
conditions.is_high_reward = isHighReward(masterMask);
conditions.is_low_reward = isLowReward(masterMask);
conditions.is_high_salience = isHighSalience(masterMask);
conditions.is_low_salience = isLowSalience(masterMask);
conditions.is_face_target = isFaceTarget(masterMask);
conditions.is_nonface_target = isNonfaceTarget(masterMask);
conditions.is_image_target = is_image_target(masterMask);
conditions.is_high_probability = isHighProbability(masterMask);
conditions.is_low_probability = isLowProbability(masterMask);

end