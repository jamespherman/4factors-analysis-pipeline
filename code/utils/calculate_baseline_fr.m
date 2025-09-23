function baseline_frs = calculate_baseline_fr(session_data, event_name, window, trial_mask)
% CALCULATE_BASELINE_FR - Computes mean firing rates for a given set of trials.
%
% This function calculates the mean firing rate for each neuron over a
% specified time window, relative to a given event, for a specific subset of
% trials.
%
% INPUTS:
%   session_data - A struct containing session-specific data, conforming to
%                  the structure defined in the project's data dictionary.
%                  It must include 'spikes' and 'eventTimes' fields.
%   event_name   - A char vector specifying the event to align to (e.g.,
%                  'fixOn'). This must be a field in session_data.eventTimes.
%   window       - A [1x2] double array defining the time window relative to
%                  the event [start, end] in seconds.
%   trial_mask   - A logical vector where 'true' indicates which trials to
%                  include in the calculation. The length must match the
%                  number of trials in the session.
%
% OUTPUT:
%   baseline_frs - A vector (nNeurons x 1) containing the average firing
%                  rate for each neuron in spikes/sec.

% --- Input Validation ---
if ~isfield(session_data.eventTimes, event_name)
    error('calculate_baseline_fr:eventNameNotFound', ...
        'The specified event_name "%s" was not found in session_data.eventTimes.', event_name);
end
if ~islogical(trial_mask)
    error('calculate_baseline_fr:invalidMask', 'trial_mask must be a logical vector.');
end
if length(trial_mask) ~= length(session_data.eventTimes.(event_name))
    error('calculate_baseline_fr:maskLengthMismatch', ...
        'The length of trial_mask must match the number of trials.');
end

% --- Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nNeurons = height(cluster_info); % Use height on table for robustness
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;
event_times = session_data.eventTimes.(event_name);

% --- Trial Selection ---
selected_trial_indices = find(trial_mask);

if isempty(selected_trial_indices)
    fprintf('Warning in calculate_baseline_fr: No trials selected by trial_mask.\n');
    baseline_frs = zeros(nNeurons, 1);
    return;
end

% --- Firing Rate Calculation ---
baseline_frs = zeros(nNeurons, 1);
window_duration = window(2) - window(1);
total_duration = length(selected_trial_indices) * window_duration;

if total_duration <= 0
    % This can happen if the window is empty or has negative duration
    baseline_frs = zeros(nNeurons, 1);
    return;
end

for i_neuron = 1:nNeurons
    neuron_id = cluster_ids(i_neuron);
    neuron_spike_times = all_spike_times(all_spike_clusters == neuron_id);

    if isempty(neuron_spike_times)
        continue; % Firing rate is already 0
    end

    total_spike_count = 0;

    for i_trial_idx = 1:length(selected_trial_indices)
        trial_idx = selected_trial_indices(i_trial_idx);

        % Get the alignment time for the current trial
        alignment_time = event_times(trial_idx);

        % Check for valid alignment time
        if isnan(alignment_time) || isinf(alignment_time)
            continue;
        end

        % Define the absolute time window for this trial
        start_time = alignment_time + window(1);
        end_time = alignment_time + window(2);

        % Count spikes within the window for this trial
        spike_count = nnz(neuron_spike_times >= start_time & neuron_spike_times < end_time);
        total_spike_count = total_spike_count + spike_count;
    end

    % Calculate the average firing rate for this neuron
    baseline_frs(i_neuron) = total_spike_count / total_duration;
end

end
