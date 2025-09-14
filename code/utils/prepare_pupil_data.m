%% prepare_pupil_data.m
%
%   Generates preprocessed, event-aligned pupil data traces for multiple
%   specified behavioral events.
%
% INPUTS:
%   session_data         - The main data structure.
%   tokens_trial_indices - Indices for trials of the 'tokens' task.
%   alignment_events     - Cell array of event names to align to.
%
% OUTPUT:
%   aligned_pupil        - A struct with aligned pupil data, with fields
%                          for each alignment event.
%
% Author: Jules
% Date: 2025-09-08

function aligned_pupil = prepare_pupil_data(session_data, ...
    tokens_trial_indices, alignment_events)

%% Define Preprocessing Parameters
sample_rate = 100; % Hz, from display refresh rate
baseline_window = [-0.5, -0.1]; % s, relative to alignment event
deblink_threshold = -9.5; % Pupil values below this are artifacts
smoothing_window_ms = 0; % ms
smoothing_window_samples = round(smoothing_window_ms / 1000 * sample_rate);
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

    % Create time vector and initialize storage based on window size
    n_samples = round((time_window(2) - time_window(1)) * sample_rate);
    time_vector = linspace(time_window(1), time_window(2), n_samples);
    aligned_traces = nan(n_tokens_trials, n_samples);

    % Get alignment times for the current event
    % Note: For 'reward', this aligns to the *first* reward pulse
    if strcmp(event_name, 'reward')
        event_times = cellfun(@(c) c(1), ...
            session_data.eventTimes.rewardCell(tokens_trial_indices));
    else
        event_times = session_data.eventTimes.(event_name)( ...
            tokens_trial_indices);
    end

    % Loop through each trial to preprocess and align pupil data
    for i_trial = 1:n_tokens_trials
        trial_idx = tokens_trial_indices(i_trial);

        if isempty(session_data.trialInfo.eyeP{trial_idx}) || ...
           isempty(session_data.trialInfo.eyeT{trial_idx})
            continue;
        end

        % pupil data 'time' is on the VIEWPixx clock so we subtract 
        % pdsTrialStartDP and add 'trialBegin' before subtracting the 
        % 'event_times{i_trial}'.
        pupil_trace = session_data.trialInfo.eyeP{trial_idx};
        pupil_time = session_data.trialInfo.eyeT{trial_idx} - ...
            session_data.eventTimes.pdsTrialStartDP(trial_idx) + ...
            session_data.eventTimes.trialBegin(trial_idx) - ...
            event_times(i_trial);

        pupil_trace(pupil_trace < deblink_threshold) = nan;

        % No smoothing if window is 0 or 1
        if smoothing_window_samples > 1
            smoothed_trace = movmean(pupil_trace, ...
                smoothing_window_samples, 'omitnan');
        else
            smoothed_trace = pupil_trace; 
        end

        baseline_indices = pupil_time >= baseline_window(1) & ...
            pupil_time <= baseline_window(2);
        baseline_mean = mean(smoothed_trace(baseline_indices), 'omitnan');

        if isnan(baseline_mean) || baseline_mean == 0
            continue;
        end

        normalized_trace = (smoothed_trace - baseline_mean) / ...
            baseline_mean;

        % Interpolate to the common time base for this event
        aligned_traces(i_trial, :) = interp1(pupil_time, ...
            normalized_trace, time_vector, 'linear', nan);
    end

    % --- Special Handling for Variable-Length 'reward' Epoch ---
    if strcmp(event_name, 'reward')
        last_reward_times_abs = cellfun(@(c) c(end), ...
            session_data.eventTimes.rewardCell(tokens_trial_indices), ...
            'UniformOutput', false);

        empty_trials = cellfun('isempty', last_reward_times_abs);
        last_reward_times_abs(empty_trials) = {NaN};
        last_reward_times_abs = cell2mat(last_reward_times_abs);

        for i_trial = 1:n_tokens_trials
            if isnan(event_times(i_trial)) || isnan(...
                    last_reward_times_abs(i_trial))
                continue;
            end

            last_reward_rel_time = last_reward_times_abs(i_trial) - ...
                event_times(i_trial);
            nan_cutoff_time = last_reward_rel_time + 1.0;
            bins_to_nan = time_vector > nan_cutoff_time;
            aligned_traces(i_trial, bins_to_nan) = nan;
        end
    end

    % Store the final matrix and its time vector
    aligned_pupil.(event_name).traces = aligned_traces;
    aligned_pupil.(event_name).time_vector = time_vector;
    aligned_pupil.(event_name).window = time_window;
end

end
