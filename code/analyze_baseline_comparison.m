%% analyze_baseline_comparison.m
%
% Compares post-event firing rates to a pre-event baseline period.
%
% This function is a plan-driven component of the 4factors analysis pipeline.
% It is designed to be called by a higher-level script that specifies the
% exact conditions and alignment events to be analyzed. For each specified
% alignment event, the function performs a statistical comparison between a
% pre-event baseline and a post-event window for a single condition.
%
% Inputs:
%   core_data        - A structure containing the core neuronal and behavioral
%                      data, pre-aligned to various events. See the project's
%                      data dictionary for more details.
%   conditions       - A structure containing logical masks for different
%                      trial conditions. The field names of this struct
%                      correspond to condition names.
%   baseline_plan_item - A structure defining the specific baseline condition
%                        to be analyzed. It must contain a 'name' field, which
%                        specifies the condition name (e.g., 'is_reward_high').
%                        This name must correspond to a field in the
%                        'conditions' struct.
%   alignment_events - A cell array of character vectors, where each element
%                      is the name of an alignment event to be processed
%                      (e.g., {'cue_onset', 'outcome_onset'}). These events
%                      must exist as fields in core_data.spikes.
%
% Outputs:
%   analysis_results - A nested structure containing the analysis results.
%                      The structure has the following hierarchy:
%                      analysis_results.EVENT_NAME.CONDITION_NAME.FIELD
%                      - EVENT_NAME: A field for each alignment event processed.
%                      - CONDITION_NAME: The name of the condition analyzed.
%                      - FIELD:
%                        - sig: A matrix (neurons x time bins) of p-values
%                          from the statistical comparison.
%                        - time_vector: The time vector for the bins in 'sig'.
%
% Author: Jules
% Date: 2025-09-20

function analysis_results = analyze_baseline_comparison(core_data, conditions, baseline_plan_item, alignment_events, stability_info)
%% Setup
% Add the 'utils' directory to the path
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

% Handle optional stability_info argument for neuron present trials
if nargin < 5
    stability_info = [];
end

% Get condition name from the baseline plan item
condition_name = baseline_plan_item.name;

% Initialize the output structure
analysis_results = struct();

%% Main Analysis Loop
% Iterate over each alignment event specified in the plan
for i_event = 1:numel(alignment_events)
    event_name = alignment_events{i_event};

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

    %% Compute neuron present trial masks if stability_info is provided
    use_stable_trials = false;
    neuron_trial_masks = [];
    if ~isempty(stability_info) && isfield(core_data, 'trial_indices')
        % First get masks for all core_data trials
        all_trial_masks = get_neuron_present_trials(stability_info, ...
            core_data.trial_indices, n_neurons);
        % Then filter to just the trials selected by trial_mask
        neuron_trial_masks = all_trial_masks(:, trial_mask);
        use_stable_trials = true;
    end

    % Initialize results arrays
    n_comp_bins = sum(comp_idx);
    sig_vals = NaN(n_neurons, n_comp_bins);

    % Iterate through each neuron
    for i_neuron = 1:n_neurons
        % Apply neuron-specific stable trial mask (if available)
        if use_stable_trials && ~isempty(neuron_trial_masks)
            stable_mask = neuron_trial_masks(i_neuron, :);
            neuron_data = sub_data(i_neuron, stable_mask, :);
            n_neuron_trials = sum(stable_mask);
        else
            neuron_data = sub_data(i_neuron, :, :);
            n_neuron_trials = n_sel_trials;
        end

        % Get baseline activity for this neuron. Squeeze to get a
        % (trials x bins) matrix:
        neuron_baseline_data = squeeze(neuron_data(1, :, baseline_idx));
        neuron_baseline_data = repmat(neuron_baseline_data(:), 1, ...
            n_comp_bins);

        % Get comparison activity for this neuron
        neuron_comp_data = squeeze(neuron_data(1, :, comp_idx));

        if isempty(neuron_baseline_data) || ...
                isempty(neuron_comp_data) || n_neuron_trials < 5
            continue;
        end

        % Compare baseline to each comparison bin
        try
            [~, ~, ~, sig_vals(i_neuron, :)] = ...
                arrayROC(neuron_baseline_data, ...
                neuron_comp_data, 500);
        catch ME
            % This can fail if one input has no variance.
            % Silently continue, leaving NaN.
        end
    end

    % Store results in the output structure
    analysis_results.(event_name).(condition_name).sig = sig_vals;
    analysis_results.(event_name).(condition_name).time_vector = ...
        bin_centers(comp_idx);
end
end
