function [conditions, condition_defs] = define_task_conditions(varargin)
% DEFINE_TASK_CONDITIONS Creates trial condition masks and defines
% analysis plans for a two-model approach.
%
% This function is the single source of truth for analysis
% configuration in the 4factors task. It partitions analyses into
% two models to de-confound identity and salience factors:
% 1.  **Image Trials Model:** Analyzes trials where the target was a
%     face or non-face image. Factors: reward, probability, identity.
% 2.  **Bullseye Trials Model:** Analyzes trials with a bullseye
%     target. Factors: reward, probability, salience.
%
% Spatial probability is determined by a fixed, block-contingent
% experimental design (e.g., Loc 1 is high-prob in Block 1), not by a
% data-driven calculation from session behavior.
%
% It returns a comprehensive analysis plan (`condition_defs`) and,
% when provided with session data, a struct of logical trial masks
% (`conditions`) aligned to this two-model approach.
%
% The `anova_plan` and `roc_plan` are now event-agnostic. The analyses
% they define will be run for all events in `condition_defs.events`.
%
% OUTPUTS:
%   conditions:     A struct of logical masks for trial conditions,
%                   specific to a session's data. It includes masks
%                   like `is_image_target` and `is_bullseye_target`
%                   to separate trials for the two models. It also
%                   contains a `factors` sub-struct with categorical
%                   labels for use in ANOVA, and a `dvs` sub-struct
%                   with pre-calculated dependent variables.
%
%   condition_defs: A struct containing the session-agnostic
%                   analysis plan. This plan is now structured as
%                   a two-element array for ANOVA and behavioral
%                   analyses, corresponding to the Image and
%                   Bullseye models.

%% --- I. Define the Comprehensive Analysis Plan ---
% This section defines the entire analysis plan. It is used by the
% main pipeline script, `run_4factors_analysis.m`, to determine
% which analyses to run.

% A. Events and Field Names
% The 'events' field is the master list of alignment events for all
% time-resolved analyses.
condition_defs.events = {'fixOn', 'targetOn', 'fixOff', 'saccadeOnset', 'reward'};

% The 'event_field_names' field specifies which struct fields in the
% analysis plans contain event names. This allows for centralized
% control over event discovery.
condition_defs.event_field_names = {'event', 'train_event', 'test_event'};

% Neuron Inclusion Configuration
% Controls whether strict screening criteria are applied to select neurons.
% When use_strict_screening is false, all kilosort clusters are included
% in analyses, which is useful for exploratory analysis.
% When true, neurons must pass ALL THREE criteria:
%   1. Mean firing rate threshold (brain-area specific, computed across entire session)
%   2. Task modulation (significant deviation from baseline at ANY location/event)
%   3. Sparsity threshold (proportion of empty 100ms bins must be < 0.7)
condition_defs.neuron_inclusion = struct(...
    'use_strict_screening', true, ...    % true = apply criteria; false = include all
    'min_firing_rate_sc', 3.5, ...       % sp/s mean FR threshold for SC neurons (changed from 5)
    'min_firing_rate_snc', 1, ...        % sp/s mean FR threshold for SNc neurons
    'max_proportion_empty_bins', 0.7, ...% max fraction of 100ms bins with 0 spikes
    'consecutive_bins_required', 3 ...   % bins needed for modulation detection (changed from 2)
);

% Modulation Windows for Task Modulation Detection
% Each event specifies baseline and test windows (in seconds) for detecting
% task-related modulation. A neuron is "modulated" if N consecutive bins
% (see consecutive_bins_required above) show significant deviation from baseline.
% For SNc neurons, all five epochs are checked: fixOn, targetOn, fixOff, saccadeOnset, reward.
condition_defs.modulation_windows = struct(...
    'fixOn', struct('baseline', [-0.25, 0], 'test', [0, 0.5]), ...
    'targetOn', struct('baseline', [-0.25, 0], 'test', [0, 0.5]), ...
    'fixOff', struct('baseline', [-0.25, 0], 'test', [0, 0.375]), ...
    'saccadeOnset', struct('baseline', [-0.35, -0.25], 'test', [-0.25, 0.25]), ...
    'reward', struct('baseline', [-0.25, 0], 'test', [0, 0.5]) ...
);

% Baseline Firing Rate Configuration
% Defines the parameters for calculating the baseline firing rate, which is
% used as a standard reference for various analyses.
condition_defs.baseline_fr_config = struct(...
    'task_name', 'gSac_4factors', ...
    'event',     'fixOn', ...
    'window',    [-0.25, 0]);

% Defines the canonical analysis window for each event
condition_defs.event_windows = struct(...
    'fixOn',        [-0.25 0.5], ...
    'targetOn',     [-0.25 0.5], ...
    'fixOff',       [-0.375, 0.375], ...
    'saccadeOnset', [-0.375, 0.375], ...
    'reward',       [-0.25 0.5] ...
);

% B. Canonical Names for Condition Masks
% These are the building blocks for the analysis plans below.
% IMPORTANT: Order is {cond1, cond2} where cond1 = "low" and cond2 = "high".
% This ensures that when passed to arrayROC(cond1, cond2), AUC > 0.5 means
% higher firing for the "high" condition (since AUC = P(X2 > X1)).
condition_defs.condition_masks.reward = ...
    {'is_low_reward', 'is_high_reward'};
condition_defs.condition_masks.salience = ...
    {'is_low_salience', 'is_high_salience'};
condition_defs.condition_masks.identity = ...
    {'is_nonface_target', 'is_face_target'};
condition_defs.condition_masks.probability = ...
    {'is_low_probability', 'is_high_probability'};
condition_defs.condition_masks.spatial = {'is_contralateral_target', ...
    'is_ipsilateral_target', 'is_opposite_rf'};

% C. Baseline Comparison Plan
% Each element defines a "Baseline vs. Post-Event Activity" analysis.
% This plan has been verified to ensure it only contains conditions
% relevant to the 4-factors task.
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

% D. ROC Analysis Plan (Bin-by-Bin)
% Each element defines a bin-by-bin ROC comparison. The ROC
% comparisons will be run for all events defined in
% condition_defs.events.
roc_plan_def = { ...
    'reward', ...
    condition_defs.condition_masks.reward{1}, ...
    condition_defs.condition_masks.reward{2}, ...
    'is_contralateral_target'; ...
    'salience', ...
    condition_defs.condition_masks.salience{1}, ...
    condition_defs.condition_masks.salience{2}, ...
    'is_contralateral_target'; ...
    'identity', ...
    condition_defs.condition_masks.identity{1}, ...
    condition_defs.condition_masks.identity{2}, ...
    'is_contralateral_target'; ...
    'probability', ...
    condition_defs.condition_masks.probability{1}, ...
    condition_defs.condition_masks.probability{2}, ...
    'is_contralateral_target' ...
};
condition_defs.roc_plan = struct('name', {}, 'cond1', {}, 'cond2', {}, 'trial_mask', {});
for i = 1:size(roc_plan_def, 1)
    condition_defs.roc_plan(i).name       = roc_plan_def{i, 1};
    condition_defs.roc_plan(i).cond1      = roc_plan_def{i, 2};
    condition_defs.roc_plan(i).cond2      = roc_plan_def{i, 3};
    condition_defs.roc_plan(i).trial_mask = roc_plan_def{i, 4};
end

% D2. Window-Based ROC Analysis Plan
% Computes ROC AUC values in fixed time windows for each neuron,
% providing a single summary metric per neuron × factor × epoch
% suitable for scatter plots and population-level comparisons.
%
% Epoch windows vary by brain area to account for different response
% latencies (SC has faster visual responses than SNc).

% Define epoch windows per brain area
% Format: struct with fields for each epoch, containing [start, end] in seconds
condition_defs.window_roc_plan.epoch_windows = struct();

% Visual epoch: aligned to targetOn
condition_defs.window_roc_plan.epoch_windows.visual = struct(...
    'event', 'targetOn', ...
    'SC',    [0.050, 0.250], ...  % 50-250ms for SC (faster visual responses)
    'SNc',   [0.100, 0.300] ...   % 100-300ms for SNc (slower visual responses)
);

% Delay epoch: aligned to fixOff (end of delay period)
condition_defs.window_roc_plan.epoch_windows.delay = struct(...
    'event', 'fixOff', ...
    'SC',    [-0.200, 0.000], ... % -200 to 0ms relative to fixOff
    'SNc',   [-0.200, 0.000] ...  % Same window for both areas
);

% Peri-saccadic epoch: aligned to saccadeOnset
condition_defs.window_roc_plan.epoch_windows.perisaccade = struct(...
    'event', 'saccadeOnset', ...
    'SC',    [-0.050, 0.100], ... % -50 to 100ms around saccade onset
    'SNc',   [-0.050, 0.100] ...  % Same window for both areas
);

% Post-reward epoch: aligned to reward
condition_defs.window_roc_plan.epoch_windows.postreward = struct(...
    'event', 'reward', ...
    'SC',    [0.100, 0.400], ...  % 100-400ms after reward
    'SNc',   [0.100, 0.400] ...   % Same window for both areas
);

% Define factor comparisons for window ROC
% Each factor specifies:
%   - name: factor identifier (reuses roc_plan names)
%   - cond1, cond2: condition mask names for high vs low comparison
%   - trial_mask: cell array of masks to AND together for trial selection
%                 This allows combining spatial restriction with factor-specific subsets
condition_defs.window_roc_plan.factors = struct(...
    'name', {}, 'cond1', {}, 'cond2', {}, 'trial_mask', {} ...
);

% Reward: all contralateral trials
condition_defs.window_roc_plan.factors(1) = struct(...
    'name', 'reward', ...
    'cond1', condition_defs.condition_masks.reward{1}, ...   % is_low_reward
    'cond2', condition_defs.condition_masks.reward{2}, ...   % is_high_reward
    'trial_mask', {{'is_contralateral_target'}} ...          % All contralateral trials
);

% Salience: bullseye trials only (contralateral)
condition_defs.window_roc_plan.factors(2) = struct(...
    'name', 'salience', ...
    'cond1', condition_defs.condition_masks.salience{1}, ... % is_low_salience
    'cond2', condition_defs.condition_masks.salience{2}, ... % is_high_salience
    'trial_mask', {{'is_contralateral_target', 'is_bullseye_target'}} ...
);

% Probability: all contralateral trials
condition_defs.window_roc_plan.factors(3) = struct(...
    'name', 'probability', ...
    'cond1', condition_defs.condition_masks.probability{1}, ... % is_low_probability
    'cond2', condition_defs.condition_masks.probability{2}, ... % is_high_probability
    'trial_mask', {{'is_contralateral_target'}} ...
);

% Identity: image trials only (contralateral)
condition_defs.window_roc_plan.factors(4) = struct(...
    'name', 'identity', ...
    'cond1', condition_defs.condition_masks.identity{1}, ... % is_nonface_target
    'cond2', condition_defs.condition_masks.identity{2}, ... % is_face_target
    'trial_mask', {{'is_contralateral_target', 'is_image_target'}} ...
);

% Analysis parameters
condition_defs.window_roc_plan.params = struct(...
    'min_trials', 5, ...           % Minimum trials per condition
    'n_bootstrap', 200, ...        % Bootstrap replicates for CI
    'alpha', 0.05 ...              % Significance level for CI
);

% Location selection mode for window ROC analysis
% 'preferred' - Use per-neuron preferred location (SC only, maximizes sensitivity)
% 'contralateral' - Pool all contralateral locations (legacy behavior)
% When 'preferred' is used:
%   - For reward/salience/identity: uses each neuron's RF-optimal location (1-4)
%   - For probability: uses location 1 or 3 (whichever has stronger visual response)
%   - is_contralateral_target is NOT applied (location handled per-neuron)
condition_defs.window_roc_plan.location_mode = 'preferred';

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

% F. N-way ANOVA Plan (Two-Model Approach)
% The ANOVA models defined here will be run for each event in
% condition_defs.events.
condition_defs.anova_plan = struct('name', {}, 'run', {}, ...
    'factors', {}, 'trial_mask', {});
% Model 1: Image Trials
condition_defs.anova_plan(1).name = 'anova_imagetrials';
condition_defs.anova_plan(1).run = true;
condition_defs.anova_plan(1).factors = {'reward', 'probability', 'identity'};
condition_defs.anova_plan(1).trial_mask = ...
    {'is_contralateral_target', 'is_image_target'};
% Model 2: Bullseye Trials
condition_defs.anova_plan(2).name = 'anova_bullseyetrials';
condition_defs.anova_plan(2).run = true;
condition_defs.anova_plan(2).factors = {'reward', 'probability', 'salience'};
condition_defs.anova_plan(2).trial_mask = ...
    {'is_contralateral_target', 'is_bullseye_target'};

% G. Behavioral Analysis Plan (Two-Model Approach)
condition_defs.behavior_plan = struct('name', {}, ...
    'dependent_variable', {}, 'factors', {}, 'trial_mask', {});
% -- Reaction Time --
condition_defs.behavior_plan(1).name = 'reaction_time_imagetrials';
condition_defs.behavior_plan(1).dependent_variable = ...
    {'eventTimes.pdsSaccadeOnset', '-', 'eventTimes.pdsFixOff'};
condition_defs.behavior_plan(1).factors = {'reward', 'probability', 'identity'};
condition_defs.behavior_plan(1).trial_mask = ...
    {'is_contralateral_target', 'is_image_target'};
condition_defs.behavior_plan(2).name = 'reaction_time_bullseyetrials';
condition_defs.behavior_plan(2).dependent_variable = ...
    {'eventTimes.pdsSaccadeOnset', '-', 'eventTimes.pdsFixOff'};
condition_defs.behavior_plan(2).factors = {'reward', 'probability', 'salience'};
condition_defs.behavior_plan(2).trial_mask = ...
    {'is_contralateral_target', 'is_bullseye_target'};
% -- Peak Velocity --
condition_defs.behavior_plan(3).name = 'peak_velocity_imagetrials';
condition_defs.behavior_plan(3).dependent_variable = {'trialInfo.peakVel'};
condition_defs.behavior_plan(3).factors = {'reward', 'probability', 'identity'};
condition_defs.behavior_plan(3).trial_mask = ...
    {'is_contralateral_target', 'is_image_target'};
condition_defs.behavior_plan(4).name = 'peak_velocity_bullseyetrials';
condition_defs.behavior_plan(4).dependent_variable = {'trialInfo.peakVel'};
condition_defs.behavior_plan(4).factors = {'reward', 'probability', 'salience'};
condition_defs.behavior_plan(4).trial_mask = ...
    {'is_contralateral_target', 'is_bullseye_target'};
% -- Endpoint Error --
condition_defs.behavior_plan(5).name = 'endpoint_error_imagetrials';
condition_defs.behavior_plan(5).dependent_variable = {'endpoint_error'};
condition_defs.behavior_plan(5).factors = {'reward', 'probability', 'identity'};
condition_defs.behavior_plan(5).trial_mask = ...
    {'is_contralateral_target', 'is_image_target'};
condition_defs.behavior_plan(6).name = 'endpoint_error_bullseyetrials';
condition_defs.behavior_plan(6).dependent_variable = {'endpoint_error'};
condition_defs.behavior_plan(6).factors = {'reward', 'probability', 'salience'};
condition_defs.behavior_plan(6).trial_mask = ...
    {'is_contralateral_target', 'is_bullseye_target'};

% H. Population Decoding Plan
% This plan is now a two-stage process. First, a set of classifiers
% are trained based on the 'training_plan'. Second, these trained
% models are tested on new data according to the 'testing_plan'.

% H1. Training Plan
% Defines all classifiers to be trained.
training_plan = struct('model_tag', {}, 'event', {}, ...
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
        entry.event = epoch_details.event;
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
                      'test_cond2', {}, 'event', {}, ...
                      'time_window', {}, 'train_event', {}, 'test_event', {}, ...
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
            entry.event = epoch_details.event;
            entry.time_window = epoch_details.window;
            % Fields for cross-time, leave empty
            entry.train_event = '';
            entry.test_event = '';
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
    train_item_idx = find(strcmp({training_plan.model_tag}, entry.train_model_tag));
    entry.train_event = training_plan(train_item_idx).event;
    entry.test_event = epoch_defs.Saccade.event;
    entry.test_time_window = epoch_defs.Saccade.window;
    entry.test_cond1 = factor_conds{1};
    entry.test_cond2 = factor_conds{2};
    % Fields for cross-factor, leave empty
    entry.event = '';
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
    train_item_idx = find(strcmp({training_plan.model_tag}, entry.train_model_tag));
    entry.train_event = training_plan(train_item_idx).event;
    entry.test_event = epoch_defs.Visual.event;
    entry.test_time_window = epoch_defs.Visual.window;
    entry.test_cond1 = factor_conds{1};
    entry.test_cond2 = factor_conds{2};
    % Fields for cross-factor, leave empty
    entry.event = '';
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
    entry.event = '';
    entry.time_window = [];
    entry.train_event = '';
    entry.test_event = '';
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

% Define visual field masks based on standard polar coordinates
theta = trialInfo.targetTheta / 10;
is_right_visual_field = (theta >= 0 & theta < 90) | ...
                        (theta > 270 & theta < 360);
is_left_visual_field = (theta > 90 & theta < 270);

if grid_x < 0
    scSide = 'left';
    % Left SC -> Contralateral is Right Visual Field
    contra_thetas = is_right_visual_field;
elseif grid_x > 0
    scSide = 'right';
    % Right SC -> Contralateral is Left Visual Field
    contra_thetas = is_left_visual_field;
else
    % Default to right visual field if grid_x is 0 or ambiguous
    scSide = 'unknown';
    contra_thetas = is_right_visual_field;
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
isHighReward = (trialInfo.rewardDuration > 200);
isLowReward = (trialInfo.rewardDuration < 200);

% 2. Salience (0=image, 1=high, 2=low)
isHighSalience = (trialInfo.salience == 1);
isLowSalience = (trialInfo.salience == 2);

% 3. Identity (stimType: 1=Face, 2=Non-face, >2=Bullseye)
isFaceTarget = (trialInfo.stimType == 1);
isNonfaceTarget = (trialInfo.stimType == 2);
is_image_target = isFaceTarget | isNonfaceTarget;
is_bullseye_target = (trialInfo.stimType > 2);

% 4. Spatial Probability (Block-Dependent)
% Defines spatial probability based on a fixed, block-contingent design.
% This approach ensures that probability is determined by the experimental
% setup, not calculated from the data, thus avoiding circularity.
%
% The rules are:
% - High-probability: Block 1 -> Loc 1; Block 2 -> Loc 3
% - Low-probability:  All other locations within the block.
isHighProbability = (trialInfo.blockNumber == 1 & trialInfo.targetLocIdx == 1) | ...
                    (trialInfo.blockNumber == 2 & trialInfo.targetLocIdx == 3);

isLowProbability = (trialInfo.blockNumber == 1 & ismember(trialInfo.targetLocIdx, [2, 3, 4])) | ...
                   (trialInfo.blockNumber == 2 & ismember(trialInfo.targetLocIdx, [1, 2, 4]));

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
conditions.is_bullseye_target = is_bullseye_target(masterMask);
conditions.is_high_probability = isHighProbability(masterMask);
conditions.is_low_probability = isLowProbability(masterMask);

% --- Target Location Index (for per-neuron preferred location analysis) ---
% Store the targetLocIdx (1-4) for each trial, filtered by masterMask
conditions.targetLocIdx = trialInfo.targetLocIdx(masterMask);

% --- SNc Subregion Masks ---
% These masks apply uniformly to all trials based on session metadata.
% They enable testing Hypothesis 2 about SNc subregion specialization.
n_trials = sum(masterMask);
if isfield(sessionData.metadata, 'snc_subregion')
    snc_subregion = sessionData.metadata.snc_subregion;
    conditions.is_rvmSNc = repmat(strcmp(snc_subregion, 'rvmSNc'), n_trials, 1);
    conditions.is_cdlSNc = repmat(strcmp(snc_subregion, 'cdlSNc'), n_trials, 1);
else
    % Default to false if not set (e.g., SC sessions or legacy data)
    conditions.is_rvmSNc = false(n_trials, 1);
    conditions.is_cdlSNc = false(n_trials, 1);
end

%% --- Create Categorical Factors for ANOVA ---
% This section creates the final categorical grouping variables for
% the ANOVA, filtered by the masterMask.
conditions.factors = struct();

% Identity Factor
identity_factor = cell(size(masterMask));
identity_factor(trialInfo.stimType == 1) = {'face'};
identity_factor(trialInfo.stimType == 2) = {'nonface'};
identity_factor(trialInfo.stimType > 2) = {'bullseye'};
conditions.factors.identity = identity_factor(masterMask);

% Salience Factor
salience_factor = cell(size(masterMask));
salience_factor(isHighSalience) = {'high'};
salience_factor(isLowSalience) = {'low'};
% For image trials, salience is not a factor, so we can mark it
% as 'neutral' or NaN. For bullseye trials, it is a key factor.
salience_factor(is_image_target) = {'neutral'};
conditions.factors.salience = salience_factor(masterMask);

% Reward Factor
reward_factor = cell(size(masterMask));
reward_factor(isHighReward) = {'high'};
reward_factor(isLowReward) = {'low'};
conditions.factors.reward = reward_factor(masterMask);

% Probability Factor
prob_factor = cell(size(masterMask));
prob_factor(isHighProbability) = {'high'};
prob_factor(isLowProbability) = {'low'};
conditions.factors.probability = prob_factor(masterMask);


%% --- Create Dependent Variables for Behavioral Analysis ---
% This section calculates all dependent variables required by the
% behavior_plan. It ensures that each DV is calculated only once and
% is correctly filtered by the masterMask.

conditions.dvs = struct();
dvs_calculated = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:length(condition_defs.behavior_plan)
    plan_item = condition_defs.behavior_plan(i);
    dv_def = plan_item.dependent_variable;

    % Derive a clean name for the DV from the plan item's name field
    % e.g., 'reaction_time_imagetrials' -> 'reaction_time'
    name_parts = strsplit(plan_item.name, '_');
    dv_name = strjoin(name_parts(1:end-1), '_');

    % If we have already calculated this DV, skip to the next
    if isKey(dvs_calculated, dv_name)
        continue;
    end

    % Calculate the full-length DV vector
    dv_full = NaN(height(trialInfo), 1);

    if strcmp(dv_name, 'endpoint_error')
        % Special case: Endpoint Error
        post_sac_xy = get_nested_field(sessionData, 'trialInfo.postSacXY');
        targ_x = get_nested_field(sessionData, 'trialInfo.targDegX');
        targ_y = get_nested_field(sessionData, 'trialInfo.targDegY');
        min_len = min([size(post_sac_xy, 1), length(targ_x), length(targ_y)]);
        dx = post_sac_xy(1:min_len, 1) - targ_x(1:min_len);
        dy = post_sac_xy(1:min_len, 2) - targ_y(1:min_len);
        dv_full(1:min_len) = sqrt(dx.^2 + dy.^2);

    elseif length(dv_def) == 3 && strcmp(dv_def{2}, '-')
        % Case: Subtraction of two fields (e.g., Reaction Time)
        val1 = get_nested_field(sessionData, dv_def{1});
        val2 = get_nested_field(sessionData, dv_def{3});
        val1 = val1(:);
        val2 = val2(:);
        min_len = min(length(val1), length(val2));
        dv_full(1:min_len) = val1(1:min_len) - val2(1:min_len);
    else
        % Case: Direct extraction of a single field (e.g., Peak Velocity)
        dv_full = get_nested_field(sessionData, dv_def{1});
        dv_full = dv_full(:);
    end

    % Filter by masterMask and store in the output struct
    conditions.dvs.(dv_name) = dv_full(masterMask);
    dvs_calculated(dv_name) = true; % Mark as calculated
end

end

function val = get_nested_field(s, f)
    % Helper function to access nested fields in a struct using a string path
    parts = strsplit(f, '.');
    val = s; % Start with the base struct
    for i = 1:length(parts)
        val = val.(parts{i}); % Iteratively index into the struct
    end
end