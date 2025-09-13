%% analyze_anova.m
%
% This function performs a 3-factor ANOVA on the firing rate of each
% neuron for each time bin around three key task events. The goal is to
% understand how neuronal activity is modulated by different experimental
% factors.
%
% The ANOVA model tests for the main effects and pairwise interactions of:
%   1. RPE (Reward Prediction Error): A continuous variable representing
%      the difference between actual and expected reward. This tests if
%      neuron firing encodes reward surprise.
%   2. Reward Distribution: A categorical variable ('Normal' vs. 'Uniform')
%      representing the reward context of the trial block. This tests if
%      neurons are sensitive to the overall statistical environment.
%   3. Flicker Condition: A categorical variable representing the visual
%      stimulus properties. This tests for sensory or cognitive effects
%      related to the predictive cues.
%
% The analysis is run independently for each neuron and each time bin.
%
% INPUTS:
%   session_data - The main data structure for the session, which contains
%                  trial-by-trial information and metadata.
%   core_data    - A struct containing binned firing rates for analyzed
%                  neurons, aligned to different task events.
%   conditions   - A struct with logical masks for various trial conditions,
%                  used to build the factors for the ANOVA.
%
% OUTPUT:
%   session_data - The input session_data struct with a new field,
%                  'analysis.anova_results', which stores a structure
%                  containing the p-values for each model term, for each
%                  neuron, time bin, and alignment event.
%
% Author: Jules
% Date: 2025-09-12
%

function session_data = analyze_anova(session_data, core_data, conditions)

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Define Alignment Events
alignment_events = {'CUE_ON', 'outcomeOn', 'reward'};

%% Prepare Factors for ANOVA
% These factors are constant across all neurons and time bins. They are
% prepared once here to avoid redundant calculations inside the main loops.

% The following section carefully selects the exact set of trials used for
% the ANOVA. This is critical for ensuring that the firing rate data in
% `core_data` is correctly matched with the trial-by-trial conditions.

% Load task codes to identify 'tokens' task trials.
codes = initCodes();
tempCueFile = session_data.trialInfo.cueFile;
tempCueFile(cellfun(@isempty, tempCueFile)) = {' '};

% First, create a broad index (`is_tokens_trial_all`) of all rewarded
% trials belonging to the 'tokens' task. This matches the full set of
% trials that *could* have been used to generate the `conditions` struct.
is_tokens_trial_all = (session_data.trialInfo.taskCode == ...
    codes.uniqueTaskCode_tokens) & ...
    ~cellfun(@isempty, session_data.eventTimes.rewardCell) & ...
    ~cellfun(@isempty, session_data.trialInfo.cueFile);

% Next, create a stricter index (`is_tokens_trial`) that includes only
% those trials that had a valid cue image (i.e., not 'blank'). These are
% the actual trials that were included in the `core_data` structure and
% will be used in the ANOVA. The number of trials in this index defines
% the size of our ANOVA factors.
is_tokens_trial = is_tokens_trial_all & ...
    ~contains(tempCueFile, 'blank');
n_trials = sum(is_tokens_trial);

% Factor 1: RPE (Continuous Predictor)
% RPE is calculated as the actual reward on a given trial minus the
% mean reward across all trials in this task. This centers the predictor
% at zero and represents a simple form of reward surprise. It is treated
% as a continuous variable in the ANOVA.
rewardAmt = session_data.trialInfo.rewardAmt(is_tokens_trial);
rpe_predictor = rewardAmt - mean(rewardAmt);

% Factor 2: Distribution (Categorical, 2 Levels)
% This factor captures the reward context of the block. In 'Normal'
% blocks, the reward distribution is narrow, making outcomes more
% predictable. In 'Uniform' blocks, the distribution is broad, making
% outcomes less predictable.
% NOTE: The `conditions` struct was created based on ALL rewarded tokens
% trials. We must use `gDist` to select the subset of those trials that
% correspond to the non-blank cue trials being used in this analysis.
gDist = find(~contains(tempCueFile(is_tokens_trial_all), 'blank'));
dist_factor = cell(n_trials, 1);
dist_factor(conditions.is_normal_dist(gDist)) = {'Normal'};
dist_factor(conditions.is_uniform_dist(gDist)) = {'Uniform'};

% Factor 3: Flicker (Categorical, 4 Levels)
% This factor represents the state of the visual cue. The four levels
% capture two dimensions:
%   - Presence vs. Absence of the flicker stimulus.
%   - Certainty of the outcome, as predicted by the cue.
flicker_factor = cell(n_trials, 1);
flicker_factor(conditions.is_noflicker_certain(gDist)) = {'CertainAbsent'};
flicker_factor(conditions.is_flicker_certain(gDist)) = {'CertainPresent'};
flicker_factor(conditions.is_flicker_omitted(gDist)) = {'UncertainAbsent'};
flicker_factor(conditions.is_flicker_surprising(gDist)) = {'UncertainPresent'};

% Detect Task Variant: Check if this is an 'AV' session by looking for the
% presence of a field names 'isAVTrial' in 'trialInfo.
is_av_session = isfield(session_data.trialInfo, 'isAVTrial');

%% Main Analysis Loop
for i_event = 1:numel(alignment_events)
    event_name = alignment_events{i_event};

    % Get dimensions from core_data for this event
    [n_neurons, n_trials, n_bins] = size(...
        core_data.spikes.(event_name).rates);

    % Pre-allocate the results structure based on the session type. This
    % improves performance and code clarity.
    if is_av_session
        term_names = {'rpe', 'dist', 'flicker', 'rpe_dist', ...
            'rpe_flicker', 'dist_flicker'};
    else
        term_names = {'rpe', 'dist', 'rpe_dist'};
    end

    for i_term = 1:numel(term_names)
        p_field_name = ['p_' term_names{i_term}];
        results.(event_name).(p_field_name) = nan(n_neurons, n_bins);
    end

    % Loop through each neuron and each time bin to perform the ANOVA
    for i_neuron = 1:n_neurons
        for i_bin = 1:n_bins

            % Extract the firing rates for the current neuron and time bin
            % across all trials. Squeeze removes singleton dimensions.
            firing_rates = squeeze(core_data.spikes.(event_name).rates( ...
                i_neuron, gDist, i_bin));

            % Anovan requires at least two groups for each categorical
            % variable. If NaNs in the firing rate data cause a factor to
            % have only one level, anovan will error. We check for this
            % and skip the bin if the condition isn't met.
            valid_trials = ~isnan(firing_rates);
            if numel(unique(dist_factor(valid_trials))) < 2
                continue;
            end

            % Conditionally run the appropriate ANOVA model based on whether
            % the session includes the flicker manipulation.
            try
                if is_av_session
                    % For AV sessions, run the full 3-factor ANOVA.
                    if numel(unique(flicker_factor(valid_trials))) < 2
                        continue; % Skip if not enough flicker groups
                    end
                    [~, tbl, ~] = anovan(firing_rates, ...
                        {rpe_predictor, dist_factor, flicker_factor}, ...
                        'model', 'interaction', 'continuous', 1, ...
                        'varnames', {'rpe', 'dist', 'flicker'}, ...
                        'display', 'off');

                    % Store p-values for all 6 terms
                    results.(event_name).p_rpe(i_neuron, i_bin) = tbl{2, 7};
                    results.(event_name).p_dist(i_neuron, i_bin) = tbl{3, 7};
                    results.(event_name).p_flicker(i_neuron, i_bin) = tbl{4, 7};
                    results.(event_name).p_rpe_dist(i_neuron, i_bin) = tbl{5, 7};
                    results.(event_name).p_rpe_flicker(i_neuron, i_bin) = tbl{6, 7};
                    results.(event_name).p_dist_flicker(i_neuron, i_bin) = tbl{7, 7};

                else
                    % For main task sessions, run a simpler 2-factor ANOVA.
                    [~, tbl, ~] = anovan(firing_rates, ...
                        {rpe_predictor, dist_factor}, ...
                        'model', 'interaction', 'continuous', 1, ...
                        'varnames', {'rpe', 'dist'}, 'display', 'off');

                    % Store p-values for the 3 relevant terms
                    results.(event_name).p_rpe(i_neuron, i_bin) = tbl{2, 7};
                    results.(event_name).p_dist(i_neuron, i_bin) = tbl{3, 7};
                    results.(event_name).p_rpe_dist(i_neuron, i_bin) = tbl{4, 7};
                end

            catch ME
                % If anovan fails, the p-values for this bin remain NaN.
                % This is expected for bins with insufficient data.
                % Uncomment the line below for detailed debugging.
                % fprintf('Could not run ANOVA for neuron %d, bin %d, event %s: %s\n', ...
                %     i_neuron, i_bin, event_name, ME.message);
            end

        end % End of time bin loop
    end % End of neuron loop
end % End of event loop

%% Finalize Output
% Assign the completed results structure to a new field in the session_data
% struct. This modified struct is the final output of the function.
session_data.analysis.anova_results = results;

fprintf('N-way ANOVA analysis complete.\n');

end
