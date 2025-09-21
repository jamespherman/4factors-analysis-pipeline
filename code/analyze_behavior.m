%% analyze_behavior.m
%
% This function executes a single behavioral analysis as defined by an
% entry in the behavior_plan. It calculates the dependent variable,
% prepares a data table, fits a linear mixed-effects model, and returns
% the results.
%
% Author: Jules
% Date: 2025-09-14
%

function results = analyze_behavior(session_data, conditions, behavior_plan_item)
% ANALYZE_BEHAVIOR Fits a linear mixed-effects model to behavioral data.
%
% INPUTS:
%   session_data:       The main session data structure.
%   conditions:         A struct of logical masks for trial conditions.
%   behavior_plan_item: A single element from the analysis_plan.behavior_plan
%                       struct array, defining the analysis to be run.
%
% OUTPUTS:
%   results:            An ANOVA table from the fitted LME model.

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% --- 1. Calculate the Dependent Variable (DV) ---
dv_definition = behavior_plan_item.dependent_variable;
num_trials = height(session_data.trialInfo);
dv = NaN(num_trials, 1);

if ischar(dv_definition{1}) && strcmp(dv_definition{1}, 'endpoint_error')
    % Special case: Calculate Endpoint Error from post-saccade eye
    % position and target location.
    post_sac_xy = get_nested_field(session_data, 'trialInfo.postSacXY');
    targ_x = get_nested_field(session_data, 'trialInfo.targDegX');
    targ_y = get_nested_field(session_data, 'trialInfo.targDegY');

    min_len = min([size(post_sac_xy, 1), length(targ_x), length(targ_y)]);

    dx = post_sac_xy(1:min_len, 1) - targ_x(1:min_len);
    dy = post_sac_xy(1:min_len, 2) - targ_y(1:min_len);
    dv(1:min_len) = sqrt(dx.^2 + dy.^2);

elseif length(dv_definition) == 3 && strcmp(dv_definition{2}, '-')
    % Case: Subtraction of two fields (e.g., for Reaction Time)
    val1 = get_nested_field(session_data, dv_definition{1});
    val2 = get_nested_field(session_data, dv_definition{3});

    val1 = val1(:);
    val2 = val2(:);
    min_len = min(length(val1), length(val2));

    dv(1:min_len) = val1(1:min_len) - val2(1:min_len);
else
    % Case: Direct extraction of a single field (e.g., for Peak Velocity)
    dv = get_nested_field(session_data, dv_definition{1});
    dv = dv(:); % Ensure it's a column vector
end


%% --- 2. Prepare the Data Table for LME ---
% Get the logical trial mask for this specific analysis
trial_mask = conditions.(behavior_plan_item.trial_mask);

% Dependent Variable
tbl = table(dv(trial_mask)', 'VariableNames', {'DV'});

% Independent Variables (Factors)
factors = behavior_plan_item.factors;
try
for i = 1:length(factors)
    factor_name = factors{i};
    % Factor data is stored in trialInfo and must be categorical for LME
    factor_data = session_data.trialInfo.(factor_name);
    tbl.(factor_name) = categorical(factor_data(trial_mask));
end
catch me
    keyboard
end

% Random Effect (Session ID)
% Replicate the session_id to match the number of trials in the mask
session_id_str = session_data.metadata.unique_id;
num_valid_trials = sum(trial_mask);
session_id_col = repmat({session_id_str}, num_valid_trials, 1);
tbl.session_id = categorical(session_id_col);


%% --- 3. Fit Linear Mixed-Effects Model ---
% Dynamically construct the model formula string, e.g.,
% 'DV ~ reward * salience * identity * probability + (1|session_id)'
factor_string = strjoin(factors, ' * ');
formula = sprintf('DV ~ %s + (1|session_id)', factor_string);

% Fit the model. Suppress non-integer response warning, as reaction
% times and other measures are often floating-point values.
warning('off', 'stats:fitlme:NonIntegerVarResp');
lme = fitlme(tbl, formula);
warning('on', 'stats:fitlme:NonIntegerVarResp');

% Extract the ANOVA table containing the statistics for each factor
results = anova(lme);

end

function val = get_nested_field(s, f)
    % Helper function to access nested fields in a struct using a string path
    parts = strsplit(f, '.');
    val = s; % Start with the base struct
    for i = 1:length(parts)
        val = val.(parts{i}); % Iteratively index into the struct
    end
end