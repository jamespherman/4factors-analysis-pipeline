function [conditions, is_av_session, condition_defs] = define_task_conditions(trialInfo, eventTimes, unique_id)
% DEFINE_TASK_CONDITIONS Creates trial condition masks and defines analysis plans
%
% This function serves a dual purpose:
% 1.  When called with session data (trialInfo, eventTimes, unique_id), it
%     calculates and returns a struct of logical masks for various trial
%     conditions based on the session's data.
% 2.  When called without arguments, it returns a comprehensive analysis
%     plan in the `condition_defs` struct. This plan is the single source
%     of truth for all analyses run in the main pipeline.
%
% OUTPUTS:
%   conditions:     Struct of logical masks for trial conditions.
%   is_av_session:  Boolean, true if the session contains AV trials.
%   condition_defs: A struct containing both the canonical names for condition
%                   masks and the full, structured analysis plan.

%% --- I. Define the Comprehensive Analysis Plan ---
% This section defines the entire analysis plan. It is used by the main
% pipeline script, `run_tokens_analysis.m`, to determine which analyses
% to run.

% A. Canonical Names for Condition Masks
% These are the building blocks for the analysis plans below.
condition_defs.condition_masks.distribution = {'is_normal_dist', 'is_uniform_dist'};
condition_defs.condition_masks.RPE_normal = {'is_norm_rare_low', 'is_norm_common', 'is_norm_rare_high'};
condition_defs.condition_masks.SPE = {'is_flicker_certain', 'is_flicker_surprising', 'is_flicker_omitted', 'is_noflicker_certain'};
condition_defs.condition_masks.RPE_comparison_pair = {'is_norm_common', 'is_norm_rare_high'};
condition_defs.condition_masks.SPE_comparison_pair = {'is_flicker_certain', 'is_flicker_surprising'};

% B. ROC Analysis Plan
% Each element defines a bin-by-bin ROC comparison.
%   .name:        Unique name for the analysis (used for saving results).
%   .event:       Alignment event (e.g., 'CUE_ON', 'outcomeOn').
%   .cond1:       Name of the first condition mask.
%   .cond2:       Name of the second condition mask.
%   .is_av_only:  Boolean, true if analysis is specific to AV sessions.
roc_plan_def = { ...
    'Dist_at_Cue',    'CUE_ON',    condition_defs.condition_masks.distribution{1},      condition_defs.condition_masks.distribution{2},      false; ...
    'RPE_at_Outcome', 'outcomeOn', condition_defs.condition_masks.RPE_comparison_pair{1}, condition_defs.condition_masks.RPE_comparison_pair{2}, false; ...
    'RPE_at_Reward',  'reward',    condition_defs.condition_masks.RPE_comparison_pair{1}, condition_defs.condition_masks.RPE_comparison_pair{2}, false; ...
    'SPE_at_Outcome', 'outcomeOn', condition_defs.condition_masks.SPE_comparison_pair{1}, condition_defs.condition_masks.SPE_comparison_pair{2}, true ...
};
condition_defs.roc_plan = struct('name', {}, 'event', {}, 'cond1', {}, 'cond2', {}, 'is_av_only', {});
for i = 1:size(roc_plan_def, 1)
    condition_defs.roc_plan(i).name       = roc_plan_def{i, 1};
    condition_defs.roc_plan(i).event      = roc_plan_def{i, 2};
    condition_defs.roc_plan(i).cond1      = roc_plan_def{i, 3};
    condition_defs.roc_plan(i).cond2      = roc_plan_def{i, 4};
    condition_defs.roc_plan(i).is_av_only = roc_plan_def{i, 5};
end

% C. Baseline Comparison Plan
% Each element defines a "Baseline vs. Post-Event Activity" analysis.
%   .name:        The name of the condition mask to use.
%   .is_av_only:  Boolean, true if analysis is specific to AV sessions.
baseline_conditions = [ ...
    condition_defs.condition_masks.RPE_normal, ...
    condition_defs.condition_masks.SPE ...
];
is_av_flags = [false, false, false, true, true, true, true];

condition_defs.baseline_plan = struct('name', {}, 'is_av_only', {});
for i = 1:length(baseline_conditions)
    condition_defs.baseline_plan(i).name       = baseline_conditions{i};
    condition_defs.baseline_plan(i).is_av_only = is_av_flags(i);
end

% D. N-way ANOVA Plan
% Each element defines an N-way ANOVA to be run.
%   .event:           Alignment event (e.g., 'CUE_ON', 'outcomeOn').
%   .p_value_names:   Name of the p-value fields to aggregate.
%   .run:             Boolean, true to run the analysis.
condition_defs.anova_plan = struct('event', {}, 'p_value_names', {}, 'run', {});
events_for_anova = {'CUE_ON', 'outcomeOn', 'reward'};
p_value_fields = { ...
    'p_value_reward', 'p_value_stim_id', 'p_value_interaction', ...
    'p_value_flicker', 'p_value_flicker_x_reward' ...
};
for i = 1:length(events_for_anova)
    condition_defs.anova_plan(i).event = events_for_anova{i};
    condition_defs.anova_plan(i).p_value_names = p_value_fields;
    condition_defs.anova_plan(i).run = true;
end


%% --- II. Mode Dispatch ---
% If the function is called without arguments, return the analysis plan.
if nargin == 0
    conditions = [];
    is_av_session = NaN;
    return;
end


%% --- III. Session-Specific Condition Mask Calculation ---
% This section runs only when session data is provided as input. It
% calculates the logical masks for each condition based on the trial data.

% A. Setup and Trial Filtering
codes = initCodes;
is_av_session = isfield(trialInfo, 'isAVTrial');

is_tokens_trial = (trialInfo.taskCode == codes.uniqueTaskCode_tokens) & ...
    ~cellfun(@isempty, eventTimes.rewardCell) & ...
    ~cellfun(@isempty, trialInfo.cueFile);

trialInfo.cueFile = trialInfo.cueFile(is_tokens_trial);
trialInfo.dist = trialInfo.dist(is_tokens_trial);
trialInfo.rewardAmt = trialInfo.rewardAmt(is_tokens_trial);
if is_av_session
    trialInfo.isAVTrial = trialInfo.isAVTrial(is_tokens_trial);
end
eventTimes.reward = eventTimes.reward(is_tokens_trial);

% B. Foundational Conditions
conditions.is_familiar = contains(trialInfo.cueFile, 'fam');
conditions.is_novel = contains(trialInfo.cueFile, 'nov');
conditions.is_normal_dist = trialInfo.dist == 1;
conditions.is_uniform_dist = trialInfo.dist == 2;
conditions.is_rewarded = eventTimes.reward > 0;

% C. Reward Magnitude / RPE Conditions
reward_norm = trialInfo.rewardAmt(conditions.is_normal_dist);
norm_thresholds = prctile(reward_norm, [25 75]);
norm_p25 = norm_thresholds(1);
norm_p75 = norm_thresholds(2);

conditions.is_norm_rare_low = conditions.is_normal_dist & (trialInfo.rewardAmt <= norm_p25);
conditions.is_norm_common = conditions.is_normal_dist & (trialInfo.rewardAmt > norm_p25 & trialInfo.rewardAmt < norm_p75);
conditions.is_norm_rare_high = conditions.is_normal_dist & (trialInfo.rewardAmt >= norm_p75);

reward_unif = trialInfo.rewardAmt(conditions.is_uniform_dist);
unif_thresholds = prctile(reward_unif, [25 75]);
unif_p25 = unif_thresholds(1);
unif_p75 = unif_thresholds(2);

conditions.is_unif_low = conditions.is_uniform_dist & (trialInfo.rewardAmt <= unif_p25);
conditions.is_unif_mid = conditions.is_uniform_dist & (trialInfo.rewardAmt > unif_p25 & trialInfo.rewardAmt < unif_p75);
conditions.is_unif_high = conditions.is_uniform_dist & (trialInfo.rewardAmt >= unif_p75);

% D. Sensory Prediction Error (SPE) Conditions
if is_av_session
    conditions.is_flicker_certain = contains(trialInfo.cueFile, '_03.jpg') & trialInfo.isAVTrial == true;
    conditions.is_flicker_surprising = contains(trialInfo.cueFile, '_02.jpg') & trialInfo.isAVTrial == true;
    conditions.is_flicker_omitted = contains(trialInfo.cueFile, '_02.jpg') & trialInfo.isAVTrial == false;
    conditions.is_noflicker_certain = contains(trialInfo.cueFile, '_01.jpg');
else
    num_trials = length(conditions.is_familiar);
    conditions.is_flicker_certain = false(num_trials, 1);
    conditions.is_flicker_surprising = false(num_trials, 1);
    conditions.is_flicker_omitted = false(num_trials, 1);
    conditions.is_noflicker_certain = false(num_trials, 1);
end

% E. Key RPE x SPE Interaction Conditions
conditions.is_common_reward_no_spe = conditions.is_norm_common & conditions.is_noflicker_certain;
conditions.is_rare_high_reward_no_spe = conditions.is_norm_rare_high & conditions.is_noflicker_certain;
conditions.is_common_reward_with_spe = conditions.is_norm_common & conditions.is_flicker_surprising;
conditions.is_rare_high_reward_with_spe = conditions.is_norm_rare_high & conditions.is_flicker_surprising;

% F. Generate and Save Diagnostic Plot
fig = figure('Visible', 'off', 'Position', [100, 100, 800, 400]);
mySubPlot([1, 2, 1]);
histogram(reward_norm);
title('Normal Distribution Outcomes');
xlabel('Reward Amount');
ylabel('Count');
hold on;
xline(norm_p25, '--r', 'LineWidth', 2);
xline(norm_p75, '--r', 'LineWidth', 2);
hold off;

mySubPlot([1, 2, 2]);
histogram(reward_unif);
title('Uniform Distribution Outcomes');
xlabel('Reward Amount');
hold on;
xline(unif_p25, '--r', 'LineWidth', 2);
xline(unif_p75, '--r', 'LineWidth', 2);
hold off;

sgtitle(sprintf('Reward Distributions for Session: %s', unique_id), 'Interpreter', 'none');
figures_dir = '../figures';
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
file_name = fullfile(figures_dir, sprintf('%s_reward_distributions.pdf', unique_id));
pdfSave(file_name, [11 8.5], fig);
close(fig);

end
