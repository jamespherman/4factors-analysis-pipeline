function analysis_plan = define_analysis_plan()
% DEFINE_ANALYSIS_PLAN creates a centralized plan for pipeline analyses.
%
% This function acts as a "control panel," defining the initial set of
% analyses to be run. It makes it easy to see and modify the planned
% analyses without changing the core processing code.
%
% OUTPUT:
%   analysis_plan: A struct containing the analysis plan.

% --- Part 1: Baseline vs. Post-Event Activity Analysis ---
% This section defines which conditions to run the baseline comparison for.
% The main analysis script will need to create these combined conditions
% on the fly by taking the logical AND of the masks from the 'conditions' struct
% (e.g., conditions.is_high_reward & conditions.is_in_rf).

analysis_plan.baseline_comparison.conditions_to_run = {
    'is_in_rf', ...
    'is_high_reward_in_rf', ...
    'is_low_reward_in_rf', ...
    'is_high_salience_in_rf', ...
    'is_low_salience_in_rf', ...
    'is_face_target_in_rf', ...
    'is_non_face_target_in_rf'
    };

% --- Part 2: Bin-by-Bin ROC Comparison Analysis ---
% This section defines the pairwise comparisons for the ROC analysis.
% Each entry in the struct array specifies the two conditions to compare.
% As with the baseline analysis, the combined conditions (e.g., 'is_high_reward_in_rf')
% must be created by the main analysis script before these comparisons are run.

comparisons(1).name = 'reward_in_rf';
comparisons(1).event = 'targetOn';
comparisons(1).cond1 = 'is_high_reward_in_rf';
comparisons(1).cond2 = 'is_low_reward_in_rf';

comparisons(2).name = 'salience_in_rf';
comparisons(2).event = 'targetOn';
comparisons(2).cond1 = 'is_high_salience_in_rf';
comparisons(2).cond2 = 'is_low_salience_in_rf';

comparisons(3).name = 'identity_in_rf';
comparisons(3).event = 'targetOn';
comparisons(3).cond1 = 'is_face_target_in_rf';
comparisons(3).cond2 = 'is_non_face_target_in_rf';

comparisons(4).name = 'probability_in_rf';
comparisons(4).event = 'targetOn';
comparisons(4).cond1 = 'is_high_probability_in_rf';
comparisons(4).cond2 = 'is_low_probability_in_rf';

analysis_plan.roc_comparison.comparisons_to_run = comparisons;

end
