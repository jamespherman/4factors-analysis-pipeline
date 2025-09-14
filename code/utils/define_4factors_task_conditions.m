function output = define_4factors_task_conditions(varargin)
% DEFINE_TASK_CONDITIONS is a dual-purpose function for the 4factors task.
%
% This function serves two roles:
% 1. Analysis Plan Definition (nargin == 0): When called with no arguments,
%    it returns the complete analysis plan, defining alignment events,
%    comparisons, and conditions for analysis.
% 2. Condition Mask Generation (nargin > 0): When called with session_data,
%    it generates logical masks for various experimental conditions based on
%    the trial data.
%
%   Usage:
%   % To get the analysis plan:
%   analysisPlan = analysis.define_task_conditions();
%
%   % To get condition masks for a session:
%   conditions = analysis.define_task_conditions(session_data);
%
%   Inputs:
%   varargin{1} - (Optional) A struct of session_data. If not provided,
%                 the function returns the analysis plan.
%
%   Outputs:
%   output - If nargin == 0, a struct named 'analysisPlan'.
%            If nargin > 0, a struct named 'conditions'.
%
% See also: analysis.determine_rf_location

if nargin == 0
    % --- Define and return the analysis plan ---
    analysisPlan.alignment_events = {'fixOn', 'targetOn', 'fixOff', ...
        'saccadeOnset', 'reward'};

    % Define conditions for "Baseline vs. Post-Event Activity" analysis
    analysisPlan.baseline_comparison.conditions_to_run = {
        'is_high_reward', 'is_low_reward', 'is_high_salience', ...
        'is_low_salience', 'is_face_target', 'is_nonface_target', ...
        'is_high_probability', 'is_low_probability'
        };

    % Define comparisons for ROC analysis
    comparisons(1).name = 'reward';
    comparisons(1).event = 'targetOn';
    comparisons(1).cond1 = 'is_high_reward';
    comparisons(1).cond2 = 'is_low_reward';
    comparisons(1).trial_mask = 'is_in_rf';

    comparisons(2).name = 'salience';
    comparisons(2).event = 'targetOn';
    comparisons(2).cond1 = 'is_high_salience';
    comparisons(2).cond2 = 'is_low_salience';
    comparisons(2).trial_mask = 'is_in_rf';

    comparisons(3).name = 'identity';
    comparisons(3).event = 'targetOn';
    comparisons(3).cond1 = 'is_face_target';
    comparisons(3).cond2 = 'is_nonface_target';
    comparisons(3).trial_mask = 'is_in_rf';

    comparisons(4).name = 'probability';
    comparisons(4).event = 'targetOn';
    comparisons(4).cond1 = 'is_high_probability';
    comparisons(4).cond2 = 'is_low_probability';
    comparisons(4).trial_mask = 'is_in_rf';

    analysisPlan.roc_comparison.comparisons_to_run = comparisons;

    % --- Per-Neuron Diagnostic Plots ---
    % Defines the set of plots to be generated for each neuron in the
    % summary PDF. Each entry in the struct array defines one panel.
    diag_plots(1).event = 'targetOn';
    diag_plots(1).title = 'Target On';
    diag_plots(1).conditions_to_compare = {'is_in_rf', 'is_out_of_rf'};

    diag_plots(2).event = 'saccadeOnset';
    diag_plots(2).title = 'Saccade Onset';
    diag_plots(2).conditions_to_compare = {'is_in_rf', 'is_out_of_rf'};

    analysisPlan.diagnostic_plots = diag_plots;

    % --- N-way ANOVA Plan ---
    analysisPlan.anova_plan.run = true;
    analysisPlan.anova_plan.event = 'targetOn';
    analysisPlan.anova_plan.factors = {'reward', 'salience', 'identity', 'probability'};
    analysisPlan.anova_plan.trial_mask = 'is_in_rf';

    output = analysisPlan;
else
    % --- Calculate and return session-specific condition masks ---
    sessionData = varargin{1};
    trialInfo = sessionData.trialInfo;

    % --- Initial Trial Filtering (The Master Mask) ---
    % Identify trials belonging to the gSac_4factors task
    isGSac4factors = trialInfo.taskCode == ...
        codes.uniqueTaskCode_gSac_4factors & ...
        ~isempty(trialInfo.pdsReward);


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

    output = conditions;
end

end
