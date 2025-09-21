%% analyze_anova.m
%
% Performs a 4-way ANOVA on neuronal firing rates for the 4-factors task.
% This function uses a fixed 4-way interaction model with the factors:
% Saliency, Identity, Reward, and Probability.
%
% The 'saliency' and 'identity' factors are constructed as three-level
% categorical variables:
% - Saliency: 'high', 'low', or 'neutral'. 'Neutral' saliency corresponds
%   to trials where a bullseye was shown instead of an image.
% - Identity: 'face', 'nonface', or 'bullseye'. 'Bullseye' corresponds to
%   trials where no face or non-face image was presented.
%
% The 'reward' and 'probability' factors are two-level categorical
% variables ('high' vs. 'low').
%
% Before running the ANOVA, the function filters the data to ensure that
% only trials with valid, non-empty labels for all four factors are
% included. This prevents errors from incomplete trial data.
%
% Inputs:
%   session_data - The main data struct for the session.
%   core_data    - The core data structure containing aligned neural data.
%   conditions   - A struct with logical masks for trial conditions.
%   anova_plan   - A struct defining the ANOVA parameters, including the
%                  alignment event and the trial mask.
%
% Output:
%   session_data - The updated session_data struct with ANOVA results.
%
% Author: Jules
% Date: 2025-09-20
%

function session_data = analyze_anova(session_data, core_data, ...
    conditions, anova_plan)

    %% Setup Paths
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'utils'));

    %% Unpack ANOVA Parameters
    event_name = anova_plan.event;
    factors = anova_plan.factors;
    trial_mask_name = anova_plan.trial_mask;

    %% Initialize Results Structure
    if ~isfield(session_data, 'analysis')
        session_data.analysis = struct();
    end
    session_data.analysis.anova_results.(event_name) = struct();

    %% Prepare Data for ANOVA
    if ~isfield(core_data.spikes, event_name)
        warning('analyze_anova:eventNotFound', ...
            'Event ''%s'' not found in core_data. Skipping ANOVA.', ...
            event_name);
        return;
    end
    event_data = core_data.spikes.(event_name);
    time_vector = event_data.time_vector;

    % Calculate mean firing rate in the post-event window (t >= 0)
    post_event_bins = time_vector >= 0;
    mean_fr = mean(event_data.rates(:, :, post_event_bins), 3, 'omitnan');

    % Get the trial mask to select trials for the ANOVA
    if ~isfield(conditions, trial_mask_name)
        warning('analyze_anova:maskNotFound', ...
            'Trial mask ''%s'' not found in conditions. Skipping ANOVA.', ...
            trial_mask_name);
        return;
    end
    trial_mask = conditions.(trial_mask_name);
    filtered_fr = mean_fr(:, trial_mask);

    %% Prepare Factors for ANOVA
    % This section creates the grouping variables for the anovan function.
    % The 'saliency' and 'identity' factors are defined as 3-level
    % categorical variables, while 'reward' and 'probability' remain
    % as 2-level factors.

    num_trials = sum(trial_mask);
    saliency = cell(num_trials, 1);
    identity = cell(num_trials, 1);
    reward = cell(num_trials, 1);
    probability = cell(num_trials, 1);

    % Saliency: high, low, neutral
    saliency(conditions.is_high_salience(trial_mask)) = {'high'};
    saliency(conditions.is_low_salience(trial_mask)) = {'low'};
    saliency(~conditions.is_image_target(trial_mask)) = {'neutral'};

    % Identity: face, nonface, bullseye
    identity(conditions.is_face_target(trial_mask)) = {'face'};
    identity(conditions.is_nonface_target(trial_mask)) = {'nonface'};
    identity(~conditions.is_image_target(trial_mask)) = {'bullseye'};

    % Reward: high, low
    reward(conditions.is_high_reward(trial_mask)) = {'high'};
    reward(conditions.is_low_reward(trial_mask)) = {'low'};

    % Probability: high, low
    probability(conditions.is_high_probability(trial_mask)) = {'high'};
    probability(conditions.is_low_probability(trial_mask)) = {'low'};

    n_neurons = size(filtered_fr, 1);
    p_values = cell(1, n_neurons);
    term_names_collected = false;
    tbl_term_names = {};

    %% Create a master mask for valid trials
    % This ensures that only trials with complete data across all factors
    % are included in the ANOVA.
    valid_trials_mask = ~cellfun('isempty', saliency) & ...
                        ~cellfun('isempty', identity) & ...
                        ~cellfun('isempty', reward) & ...
                        ~cellfun('isempty', probability);

    % Filter all data by the valid trials mask
    saliency_filt = saliency(valid_trials_mask);
    identity_filt = identity(valid_trials_mask);
    reward_filt = reward(valid_trials_mask);
    probability_filt = probability(valid_trials_mask);
    fr_filt = filtered_fr(:, valid_trials_mask);

    % Prepare inputs for anovan
    group_vars = {saliency_filt, identity_filt, reward_filt, ...
                  probability_filt};
    factor_names = {'Saliency', 'Identity', 'Reward', 'Probability'};

    for i_neuron = 1:n_neurons
        y = fr_filt(i_neuron, :)';
        if sum(~isnan(y)) < length(y) * 0.5 || isempty(y)
            continue;
        end

        [p, tbl, ~, ~] = anovan(y, group_vars, 'model', 'interaction', ...
            'varnames', factor_names, 'display', 'off');
        p_values{i_neuron} = p;

        if ~term_names_collected && ~isempty(tbl)
            tbl_term_names = tbl(2:(end-2), 1);
            term_names_collected = true;
        end
    end

    %% Store Results
    if term_names_collected
        for i_term = 1:length(tbl_term_names)
            field_name = matlab.lang.makeValidName(tbl_term_names{i_term});
            session_data.analysis.anova_results.(event_name).(field_name) = NaN(1, n_neurons);
        end

        for i_neuron = 1:n_neurons
            if ~isempty(p_values{i_neuron})
                for i_term = 1:length(tbl_term_names)
                    field_name = matlab.lang.makeValidName(tbl_term_names{i_term});
                    p_val = p_values{i_neuron}(i_term);
                    session_data.analysis.anova_results.(event_name).(field_name)(i_neuron) = p_val;
                end
            end
        end
    end
end
