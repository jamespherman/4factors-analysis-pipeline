function baseline_frs = calculate_baseline_fr(session_data, neuron_ids)
% CALCULATE_BASELINE_FR Computes mean baseline firing rate for gSac_4factors task.
%
%   This function calculates a single, average baseline firing rate for each
%   specified neuron. The calculation is based exclusively on the baseline
%   periods of trials belonging to the 'gSac_4factors' task. The baseline
%   period for each trial is defined as the interval between the 'fixAq'
%   (fixation acquired) and 'targOn' (target onset) events.
%
%   Usage:
%   baseline_frs = calculate_baseline_fr(session_data, neuron_ids)
%
%   Inputs:
%   session_data - The main data struct for a session, containing trial
%                  information, event times, and spike data.
%   neuron_ids   - A numeric vector of cluster_ids for which to calculate
%                  the baseline firing rate.
%
%   Outputs:
%   baseline_frs - A column vector of the same length as neuron_ids,
%                  containing the mean baseline firing rate for each neuron
%                  in spikes per second.
%

% --- 1. Initialize Task Codes ---
% Get the standardized task code definitions.
codes = utils.initCodes();

% --- 2. Filter for gSac_4factors Trials ---
% Create a logical index to identify trials for the gSac_4factors task.
gSac_trials_idx = session_data.trialInfo.taskCode == ...
    codes.uniqueTrialCode_gSac_4factors;

% Check if any trials for the specified task were found.
if sum(gSac_trials_idx) == 0
    warning('No trials found for the gSac_4factors task. Returning NaNs.');
    baseline_frs = nan(length(neuron_ids), 1);
    return;
end

% --- 3. Subset Relevant Data ---
% Use the logical index to filter trial-based event time vectors.
fixAq_times = session_data.eventTimes.fixAq(gSac_trials_idx);
targOn_times = session_data.eventTimes.targOn(gSac_trials_idx);

% --- 4. Define Baseline Period & Calculate Durations ---
% The baseline period is the interval [fixAq_times, targOn_times).
% Calculate the duration of each baseline period.
baseline_durations = targOn_times - fixAq_times;

% --- 5. Iterate Through Neurons & Calculate Rates ---
% Initialize the output vector with NaNs.
baseline_frs = nan(length(neuron_ids), 1);

% Loop through each requested neuron_id.
for i = 1:length(neuron_ids)
    neuron_id = neuron_ids(i);

    % Find the index corresponding to the neuron_id in the session_data.
    % Assuming session_data.spikes.cluster_id is a list of all neuron ids.
    neuron_idx = find(session_data.spikes.cluster_id == neuron_id);

    if isempty(neuron_idx)
        warning('Neuron ID %d not found in session_data. Skipping.', neuron_id);
        continue;
    end

    % Get all spike times for the current neuron.
    spike_times = session_data.spikes.times{neuron_idx};

    % Initialize array to store firing rates for each trial for this neuron.
    trial_rates = [];

    % Loop through each trial in the subset.
    for j = 1:length(fixAq_times)
        start_time = fixAq_times(j);
        end_time = targOn_times(j);
        duration = baseline_durations(j);

        % Ensure the duration is positive to avoid division by zero.
        if duration > 0
            % Count spikes that fall within the baseline window [start, end).
            spike_count = sum(spike_times >= start_time & spike_times < end_time);

            % Calculate firing rate for this trial and store it.
            trial_rates(end+1) = spike_count / duration;
        end
    end

    % Calculate the mean of the trial-by-trial rates for the current neuron.
    if ~isempty(trial_rates)
        baseline_frs(i) = mean(trial_rates, 'omitnan');
    end
end

end
