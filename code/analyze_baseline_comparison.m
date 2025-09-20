%% analyze_baseline_comparison.m
%
% Compares post-event firing rates to a pre-event baseline period for a
% single specified condition.
%
% This function is designed to be a modular component of the 4factors
% analysis pipeline. It iterates through all available alignment events in
% core_data and, for a single specified condition, performs a statistical
% comparison between a pre-event baseline and a post-event window.
%
% The specific condition to be analyzed is passed in as a name-value pair.
%
% Inputs:
%   core_data     - A structure containing the core neuronal and behavioral
%                   data, pre-aligned to various events. See the project's
%                   data dictionary for more details.
%   conditions    - A structure containing logical masks for different
%                   trial conditions. The field names of this struct
%                   correspond to condition names.
%
% Name-Value Pair Inputs:
%   'condition'   - A char vector specifying the name of the condition to
%                   analyze (e.g., 'is_reward_high'). This condition must
%                   exist as a field in the 'conditions' struct.
%
% Outputs:
%   analysis_results - A structure where each field corresponds to an
%                      alignment event. Each of these fields is a struct
%                      whose name is the condition that was analyzed. This
%                      nested struct contains the following fields:
%                        - sig: A matrix (neurons x time bins) of p-values
%                          from the statistical comparison.
%                        - time_vector: The time vector for the bins in 'sig'.
%
% Author: Jules
% Date: 2025-09-20

function analysis_results = analyze_baseline_comparison(core_data, conditions, varargin)
    %% Setup
    % Add the 'utils' directory to the path
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'utils'));

    % --- Input Parsing ---
    p = inputParser;
    p.KeepUnmatched = true; % Allow other arguments not defined here
    addRequired(p, 'core_data', @isstruct);
    addRequired(p, 'conditions', @isstruct);
    addParameter(p, 'condition', '', @ischar);
    parse(p, core_data, conditions, varargin{:});

    condition_name = p.Results.condition;

    if isempty(condition_name)
        error('analyze_baseline_comparison:noCondition', ...
              'The ''condition'' name-value pair is required.');
    end

    % Initialize the output structure
    analysis_results = struct();

    % Get alignment event names dynamically from core_data.spikes
    align_events = fieldnames(core_data.spikes);

    %% Main Analysis Loop
    % Iterate over each alignment event
    for i_event = 1:numel(align_events)
        event_name = align_events{i_event};

        if ~isfield(core_data.spikes, event_name)
            continue;
        end

        event_data = core_data.spikes.(event_name);

        % Define baseline and comparison windows based on event type
        if isfield(event_data, 'pre_time')
            baseline_window = [event_data.pre_time, 0];
            if event_data.pre_time >= 0
                baseline_window = [-0.2, 0];
            end
        else
            baseline_window = [-0.2, 0]; % Default if not specified
        end

        if isfield(event_data, 'post_time')
            comp_window = [0, event_data.post_time];
            if event_data.post_time <= 0
                comp_window = [0, 0.5];
            end
        else
            comp_window = [0, 0.5]; % Default if not specified
        end


        % Get indices for baseline and comparison bins
        bin_centers = event_data.time_vector;
        baseline_idx = (bin_centers >= baseline_window(1)) & ...
                       (bin_centers < baseline_window(2));
        comp_idx = (bin_centers >= comp_window(1)) & ...
                   (bin_centers <= comp_window(2));

        if ~any(baseline_idx) || ~any(comp_idx)
            continue; % Skip if no bins in windows
        end

        % Get the trial mask for the specified condition
        if isfield(conditions, condition_name)
            trial_mask = conditions.(condition_name);
        else
            warning('analyze_baseline_comparison:conditionNotFound', ...
                'Condition ''%s'' not found in conditions struct. Skipping for event %s.', ...
                condition_name, event_name);
            continue;
        end

        if sum(trial_mask) < 2
            continue; % Skip if not enough trials for comparison
        end

        % Select data based on the trial mask
        % The data array dimensions are (neurons, trials, time_bins)
        sub_data = event_data.rates(:, trial_mask, :);
        [n_neurons, n_sel_trials, ~] = size(sub_data);

        % Initialize results arrays
        n_comp_bins = sum(comp_idx);
        sig_vals = NaN(n_neurons, n_comp_bins);

        % Iterate through each neuron
        for i_neuron = 1:n_neurons
            % Get baseline activity for this neuron
            % Squeeze to get a (trials x bins) matrix, then mean over bins
            neuron_baseline_data = mean(squeeze(sub_data(i_neuron, :, baseline_idx)), 2, 'omitnan');

            % Get comparison activity for this neuron
            neuron_comp_data = squeeze(sub_data(i_neuron, :, comp_idx));

            if n_sel_trials == 1 && n_comp_bins > 1
                neuron_comp_data = neuron_comp_data';
            end

            if isempty(neuron_baseline_data) || isempty(neuron_comp_data)
                continue;
            end

            % Compare baseline to each comparison bin
            for i_bin = 1:n_comp_bins
                comp_data_this_bin = neuron_comp_data(:, i_bin);
                if sum(~isnan(neuron_baseline_data)) < 2 || sum(~isnan(comp_data_this_bin)) < 2
                    continue;
                end
                try
                    [~, ~, ~, sig_val] = arrayROC(neuron_baseline_data, comp_data_this_bin);
                    sig_vals(i_neuron, i_bin) = sig_val;
                catch ME
                    % This can fail if one input has no variance.
                    % Silently continue, leaving NaN.
                end
            end
        end

        % Store results in the output structure
        analysis_results.(event_name).(condition_name).sig = sig_vals;
        analysis_results.(event_name).(condition_name).time_vector = bin_centers(comp_idx);
    end
end
