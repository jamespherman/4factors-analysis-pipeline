%% analyze_baseline_comparison.m
%
% Compares post-event firing rates to a pre-event baseline period for a
% single specified condition.
%
% Author: Jules
% Date: 2025-09-14

function analysis_results = analyze_baseline_comparison(core_data, condition_name)
    %% Setup
    % Add the 'utils' directory to the path
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'utils'));

    % Initialize the output structure
    analysis_results = struct();

    % Get alignment event names from core_data
    align_events = fieldnames(core_data.aligned_data);

    %% Main Analysis Loop
    % Iterate over each alignment event
    for i_event = 1:numel(align_events)
        event_name = align_events{i_event};
        event_data = core_data.aligned_data.(event_name);

        % Define baseline and comparison windows
        if strcmp(event_name, 'sac_saccadeOn')
            baseline_window = [-0.5, -0.4];
        else
            baseline_window = [event_data.pre_time, 0];
            if event_data.pre_time >= 0
                baseline_window = [-0.2, 0];
            end
        end
        comp_window = [0, event_data.post_time];
        if event_data.post_time <= 0
            comp_window = [0, 0.5];
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
        if isfield(core_data.conditions, condition_name)
            trial_mask = core_data.conditions.(condition_name);
        else
            warning('Condition not found: %s. Skipping for event %s.', ...
                    condition_name, event_name);
            continue;
        end

        if sum(trial_mask) < 2
            continue; % Skip if not enough trials
        end

        % Select data based on trial mask
        sub_data = event_data.data_array(trial_mask, :, :);
        [n_sel_trials, n_neurons, ~] = size(sub_data);

        % Initialize results arrays
        n_comp_bins = sum(comp_idx);
        roc_vals = NaN(n_neurons, n_comp_bins);
        sig_vals = NaN(n_neurons, n_comp_bins);

        % Iterate through each neuron
        for i_neuron = 1:n_neurons
            % Get baseline activity for this neuron
            neuron_baseline_data = mean(squeeze(sub_data(:, i_neuron, baseline_idx)), 2, 'omitnan');

            % Get comparison activity for this neuron
            neuron_comp_data = squeeze(sub_data(:, i_neuron, comp_idx));

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
                    [roc_val, ~, ~, sig_val] = arrayROC(neuron_baseline_data, comp_data_this_bin);
                    roc_vals(i_neuron, i_bin) = roc_val;
                    sig_vals(i_neuron, i_bin) = sig_val;
                catch ME
                    % Could add a warning here if needed
                end
            end
        end

        % Store results in the output structure
        analysis_results.(event_name).(condition_name).sig = sig_vals;
        analysis_results.(event_name).(condition_name).time_vector = bin_centers(comp_idx);
    end
end
