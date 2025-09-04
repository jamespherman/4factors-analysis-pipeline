function baseline_frs = calculate_baseline_fr(session_data)
% calculate_baseline_fr - Computes the average baseline firing rate for neurons.
%
% This function calculates the firing rate for each neuron during a
% predefined baseline period, averaged across all valid trials.
%
% INPUTS:
%   session_data - A struct containing session-specific data. It must include:
%                  - nClusters: The number of neurons.
%                  - spikes.times: A cell array of spike timestamps for each neuron.
%                  - trials.targetOnTime: A vector of target onset times for each trial.
%
% OUTPUT:
%   baseline_frs - A vector (nNeurons x 1) containing the average baseline
%                  firing rate for each neuron in spikes/sec.
%

% --- Setup ---
nNeurons = session_data.nClusters;
spike_times_all = session_data.spikes.times;

% Assumption: Trial event data is stored in a 'trials' sub-struct.
% This may need to be adjusted if the final data structure is different.
if isfield(session_data, 'trials') && isfield(session_data.trials, 'targetOnTime')
    target_on_times = session_data.trials.targetOnTime;
else
    % Fallback for gsac_data structure compatibility
    target_on_times = session_data.targetOnTime;
end


% --- Define Baseline Period ---
% We define the baseline period as the 100ms window immediately preceding
% the target onset on each trial.
baseline_start_offset = -0.100; % 100 ms before target onset
baseline_end_offset   = 0.000;  % at target onset
baseline_duration_sec = baseline_end_offset - baseline_start_offset; % should be 0.1

% --- Calculate Firing Rates ---
baseline_frs = zeros(nNeurons, 1);
valid_trials = ~isnan(target_on_times);
nValidTrials = sum(valid_trials);

if nValidTrials == 0
    fprintf('WARNING in calculate_baseline_fr: No valid trials found.\n');
    return; % Returns a vector of zeros
end

% Total duration across all valid trials
total_baseline_duration = nValidTrials * baseline_duration_sec;

for i_neuron = 1:nNeurons
    neuron_spike_times = spike_times_all{i_neuron};
    total_spike_count = 0;

    % Iterate through only the valid trials
    for i_trial = find(valid_trials)' % Loop through indices of valid trials

        event_time = target_on_times(i_trial);
        start_time = event_time + baseline_start_offset;
        end_time   = event_time + baseline_end_offset;

        % Count spikes within the baseline window for this trial
        spike_count_trial = sum(neuron_spike_times >= start_time & neuron_spike_times < end_time);
        total_spike_count = total_spike_count + spike_count_trial;
    end

    % Calculate average firing rate for this neuron
    if total_baseline_duration > 0
        baseline_frs(i_neuron) = total_spike_count / total_baseline_duration;
    else
        baseline_frs(i_neuron) = 0; % Should not happen if nValidTrials > 0
    end
end

end
