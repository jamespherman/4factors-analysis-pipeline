%% prepare_neuronal_data.m
%
%   Generates binned, event-aligned firing rate matrices for selected
%   neurons across multiple specified behavioral events.
%
% INPUTS:
%   session_data         - The main data structure.
%   selected_neurons     - A logical vector for neuron selection.
%   tokens_trial_indices - Indices for trials of the 'tokens' task.
%   alignment_events     - Cell array of event names to align to.
%
% OUTPUT:
%   aligned_spikes       - A struct with aligned spike data, with fields
%                          for each alignment event.
%
% Author: Jules
% Date: 2025-09-08

function aligned_spikes = prepare_neuronal_data(session_data, ...
    selected_neurons, tokens_trial_indices, alignment_events)

%% Define Alignment Parameters
bin_width = 0.2; % 200ms bin size
step_size = 0.1; % 100ms step size

% Get basic info
neuron_cluster_ids = find(selected_neurons);
n_selected_neurons = numel(neuron_cluster_ids);
n_tokens_trials = numel(tokens_trial_indices);

%% Process Each Alignment Event
for i_event = 1:numel(alignment_events)
    event_name = alignment_events{i_event};

    % Define time window based on the event type
    if strcmp(event_name, 'reward')
        time_window = [-0.5, 5.0]; % Extended window for reward epoch
    else
        time_window = [-0.5, 1.5]; % Standard window for other events
    end

    % --- New Overlapping Binning ---
    % Define bin start times for the sliding window
    bin_starts = time_window(1):step_size:(time_window(2) - bin_width);

    % The time vector represents the center of each bin
    time_vector = bin_starts + bin_width / 2;
    n_time_bins = numel(time_vector);

    % Get alignment times for the current event
    % Note: For 'reward', this aligns to the *first* reward pulse
    if strcmp(event_name, 'reward')
        event_times = cellfun(@(c) c(1), ...
            session_data.eventTimes.rewardCell(tokens_trial_indices));
    else
        event_times = session_data.eventTimes.(event_name)( ...
            tokens_trial_indices);
    end

    % Initialize storage for binned firing rates
    binned_rates = nan(n_selected_neurons, n_tokens_trials, n_time_bins);

    % Loop through each selected neuron to bin their spike times
    for i_neuron = 1:n_selected_neurons
        cluster_id = neuron_cluster_ids(i_neuron);
        spike_times = session_data.spikes.times( ...
            session_data.spikes.clusters == cluster_id);

        if isempty(spike_times)
            binned_rates(i_neuron, :, :) = 0;
            continue;
        end

        % Use alignAndBinSpikes to perform the binning
        [binned_counts, time_vector_out, ~] = alignAndBinSpikes(...
            spike_times, event_times, time_window(1), time_window(2), ...
            bin_width, step_size);

        % The output time vector should be consistent, but we'll use the
        % one from the function to be precise.
        if i_neuron == 1
             time_vector = time_vector_out;
        end

        % Convert spike counts to firing rate (spikes/sec)
        % and reshape to fit the destination matrix
        binned_rates(i_neuron, :, :) = reshape(...
            binned_counts / bin_width, ...
            [1, size(binned_counts,1), size(binned_counts,2)]);
    end

    % --- Special Handling for Variable-Length 'reward' Epoch ---
    if strcmp(event_name, 'reward')
        % Get the time of the last reward pulse for each trial
        % 'UniformOutput' is false to handle trials with no rewards
        last_reward_times_abs = cellfun(@(c) c(end), ...
            session_data.eventTimes.rewardCell(tokens_trial_indices), ...
            'UniformOutput', false);

        % Replace empty values with NaN and convert to a numeric array
        empty_trials = cellfun('isempty', last_reward_times_abs);
        last_reward_times_abs(empty_trials) = {NaN};
        last_reward_times_abs = cell2mat(last_reward_times_abs);

        % Loop through each trial to apply NaN padding post-alignment
        for i_trial = 1:n_tokens_trials
            if isnan(event_times(i_trial)) || ...
                    isnan(last_reward_times_abs(i_trial))
                continue;
            end

            % Calculate last reward time relative to the alignment point
            last_reward_rel_time = last_reward_times_abs(i_trial) - ...
                event_times(i_trial);

            % Define the cutoff point for valid data (1s after last reward)
            nan_cutoff_time = last_reward_rel_time + 1.0;

            % Find all time bins that occur after this cutoff
            bins_to_nan = time_vector > nan_cutoff_time;

            % Set the data in these bins to NaN for the current trial
            binned_rates(:, i_trial, bins_to_nan) = nan;
        end
    end

    try

    % Store the final processed data for this event
    aligned_spikes.(event_name).rates = binned_rates;
    aligned_spikes.(event_name).time_vector = time_vector;
    aligned_spikes.(event_name).window = time_window;
    
    catch me
        keyboard
    end
end

end
