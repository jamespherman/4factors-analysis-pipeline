function baseline_frs = calculate_baseline_fr(session_data)
% calculate_baseline_fr - Computes baseline firing rates for the 'tokens' task.
%
% This function is specialized for the 'tokens' task. It calculates the
% mean firing rate for each neuron during a dynamically determined baseline
% period. The baseline for each trial is the window between the trial's
% start and the first subsequent event (CUE_ON, pdsOutcomeOn, or reward).
%
% INPUTS:
%   session_data - A struct containing session-specific data, conforming to
%                  the structure defined in 'docs/preprocessing_docs/session_data_dictionary.md'.
%                  It must include:
%                  - spikes.cluster_info: Table with information about each cluster.
%                  - spikes.times: Vector of all spike times.
%                  - spikes.clusters: Vector mapping each spike to a cluster ID.
%                  - eventTimes.targetOn: Vector of target onset times for each trial.
%
% OUTPUT:
%   baseline_frs - A vector (nNeurons x 1) containing the average baseline
%                  firing rate for each neuron in spikes/sec.
%

% --- Setup and Data Extraction ---
codes = initCodes();
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nNeurons = height(cluster_info.cluster_id);

all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

% --- Task-Specific Trial Selection ---
% Identify trials belonging to the 'tokens' task.
tokens_trial_indices = find(session_data.trialInfo.taskCode == ...
    codes.uniqueTaskCode_tokens);

if isempty(tokens_trial_indices)
    fprintf("Warning in calculate_baseline_fr: No 'tokens' task " + ...
        "trials found.\n");
    baseline_frs = zeros(nNeurons, 1);
    return;
end

% --- Firing Rate Calculation ---
baseline_frs = zeros(nNeurons, 1);
total_baseline_duration_per_neuron = zeros(nNeurons, 1);
total_spike_count_per_neuron = zeros(nNeurons, 1);

% Iterate through each neuron to calculate its average firing rate.
for i_neuron = 1:nNeurons
    neuron_id = cluster_ids(i_neuron);
    neuron_spike_times = all_spike_times(all_spike_clusters == neuron_id);

    if isempty(neuron_spike_times)
        continue; % Skip neurons with no spikes.
    end

    % --- Dynamic Baseline Window Calculation ---
    for i_trial = 1:length(tokens_trial_indices)
        trial_idx = tokens_trial_indices(i_trial);
        trial_start_time = session_data.eventTimes.trialBegin(trial_idx);

        % Convert event times to trial-relative time for comparison.
        % CUE_ON and reward are on the master clock; pdsOutcomeOn is 
        % already trial-relative.
        cue_on_relative = session_data.eventTimes.CUE_ON(trial_idx) - ...
        trial_start_time;
        reward_relative = session_data.eventTimes.reward(trial_idx) - ...
            trial_start_time;
        pds_outcome_relative = session_data.eventTimes.pdsOutcomeOn(...
            trial_idx);

        % Find the earliest positive event time.
        event_times = [cue_on_relative, reward_relative, ...
            pds_outcome_relative];
        positive_event_times = event_times(event_times > 0);

        if isempty(positive_event_times)
            continue; % Skip trial if no valid positive event time 
            % is found.
        end

        baseline_end_time = min(positive_event_times);
        baseline_duration = baseline_end_time; % Starts from 0 
        % (trial_start_time)

        if baseline_duration <= 0
            continue; % Skip if baseline duration is not positive.
        end

        % Define the absolute time window for spike counting.
        start_time_abs = trial_start_time;
        end_time_abs = trial_start_time + baseline_duration;

        % Count spikes within this trial's baseline window.
        spike_count = nnz(neuron_spike_times >= start_time_abs & ...
            neuron_spike_times < end_time_abs);

        total_spike_count_per_neuron(i_neuron) = total_spike_count_per_neuron(i_neuron) + spike_count;
        total_baseline_duration_per_neuron(i_neuron) = total_baseline_duration_per_neuron(i_neuron) + baseline_duration;
    end

    % --- Calculate Average Firing Rate ---
    % Avoid division by zero if a neuron had no valid baseline windows.
    if total_baseline_duration_per_neuron(i_neuron) > 0
        baseline_frs(i_neuron) = total_spike_count_per_neuron(i_neuron) / total_baseline_duration_per_neuron(i_neuron);
    end
end

end
