function neuron_trial_masks = get_neuron_present_trials(stability_info, ...
    fourfactors_trial_indices, n_neurons)
% GET_NEURON_PRESENT_TRIALS Returns per-neuron trial masks for stable recording periods.
%
%   neuron_trial_masks = GET_NEURON_PRESENT_TRIALS(stability_info, ...
%       fourfactors_trial_indices, n_neurons)
%
%   This function maps the stable trial periods identified by
%   detect_stable_trial_period.m to the filtered gSac_4factors trials used
%   in core_data. For each neuron, it returns a logical mask indicating
%   which core_data trials fall within that neuron's stable recording period.
%
%   INPUTS:
%       stability_info           - Struct array from detect_stable_trial_period
%                                  or session_data.stability_info. Must have fields:
%                                  .stable_trial_start, .stable_trial_end, .is_stable
%       fourfactors_trial_indices - [n_trials x 1] Vector of trial indices used in
%                                   core_data (indices into session_data.eventTimes)
%       n_neurons                - Number of neurons (for validation)
%
%   OUTPUT:
%       neuron_trial_masks       - [n_neurons x n_trials] logical matrix
%                                  True where trial is within neuron's stable period
%                                  If stability_info is empty or neuron has no stable
%                                  period, all trials are marked as valid (true)
%
%   NOTES:
%       - Neurons without a stable period (is_stable = false) will have all
%         trials marked as valid to maintain backward compatibility
%       - The stable_trial_start/end indices refer to ALL session trials,
%         while fourfactors_trial_indices are the filtered task trials
%
%   Example:
%       % Get trial masks for use in analysis
%       masks = get_neuron_present_trials(session_data.stability_info, ...
%           core_data.trial_indices, n_neurons);
%
%       % Use in per-neuron analysis loop
%       for i_neuron = 1:n_neurons
%           stable_trials = masks(i_neuron, :);
%           neuron_data = rates(i_neuron, stable_trials, :);
%       end
%
%   See also: detect_stable_trial_period, apply_neuron_screening
%
%   Author: Claude Code
%   Date: 2026-01-24

%% Input Validation
n_trials = length(fourfactors_trial_indices);

% Initialize output: default to all trials valid (backward compatible)
neuron_trial_masks = true(n_neurons, n_trials);

% Handle empty or missing stability_info
if isempty(stability_info)
    warning('get_neuron_present_trials:noStabilityInfo', ...
        'No stability_info provided. Returning all trials as valid.');
    return;
end

% Validate number of neurons matches
n_stability_neurons = length(stability_info);
if n_stability_neurons ~= n_neurons
    warning('get_neuron_present_trials:neuronMismatch', ...
        'Number of neurons (%d) does not match stability_info length (%d). Using min.', ...
        n_neurons, n_stability_neurons);
    n_neurons = min(n_neurons, n_stability_neurons);
end

%% Map Stable Periods to Filtered Trial Indices
for i_neuron = 1:n_neurons
    % Check if this neuron has a valid stable period
    if ~isfield(stability_info(i_neuron), 'is_stable') || ...
       ~stability_info(i_neuron).is_stable
        % No stable period - keep all trials valid (default)
        continue;
    end

    % Get stable trial range (indices into ALL session trials)
    stable_start = stability_info(i_neuron).stable_trial_start;
    stable_end = stability_info(i_neuron).stable_trial_end;

    % Handle missing or invalid values
    if isnan(stable_start) || isnan(stable_end) || ...
       stable_start < 1 || stable_end < stable_start
        continue;  % Keep all trials valid
    end

    % Create mask: which fourfactors trials fall within stable period?
    % fourfactors_trial_indices(j) is the session-level trial index for
    % the j-th trial in core_data
    stable_mask = (fourfactors_trial_indices >= stable_start) & ...
                  (fourfactors_trial_indices <= stable_end);

    neuron_trial_masks(i_neuron, :) = stable_mask;
end

%% Report Summary Statistics
n_with_stable = sum([stability_info(1:n_neurons).is_stable]);
if n_with_stable > 0
    % Calculate mean proportion of trials included per neuron with stable period
    stable_neurons = find([stability_info(1:n_neurons).is_stable]);
    mean_prop = mean(sum(neuron_trial_masks(stable_neurons, :), 2) / n_trials);
    fprintf('get_neuron_present_trials: %d/%d neurons have stable periods (mean %.1f%% trials included)\n', ...
        n_with_stable, n_neurons, mean_prop * 100);
end

end
