%% analyze_roc_comparison.m
%
% This function performs an ROC analysis to compare neural firing rates
% between two specified conditions for a single comparison event.
%
% Author: Jules
% Date: 2025-09-14

function analysis_results = analyze_roc_comparison(comparison, core_data, conditions)
    %% Setup Paths
    % Add the 'utils' directory to the path so that helper functions can be
    % found.
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'utils'));

    %% Initialize Output
    analysis_results = struct();

    %% Extract Comparison Parameters
    event_name = comparison.event;
    cond1_name = comparison.cond1;
    cond2_name = comparison.cond2;
    comp_name = comparison.name;

    % Check if the alignment event exists in the core data
    if ~isfield(core_data.aligned_data, event_name)
        warning('Event ''%s'' not found in core_data. Skipping comparison ''%s''.', ...
                event_name, comp_name);
        return;
    end

    event_data = core_data.aligned_data.(event_name);
    [~, n_neurons, n_bins] = size(event_data.data_array);

    %% Get Trial Masks
    if ~isfield(conditions, cond1_name) || ~isfield(conditions, cond2_name)
        warning('One or both conditions (''%s'', ''%s'') not found. Skipping comparison ''%s''.', ...
                cond1_name, cond2_name, comp_name);
        return;
    end

    cond1_mask = conditions.(cond1_name);
    cond2_mask = conditions.(cond2_name);

    % Apply an overall trial mask if specified
    if isfield(comparison, 'trial_mask') && ~isempty(comparison.trial_mask)
        if isfield(conditions, comparison.trial_mask)
            trial_mask = conditions.(comparison.trial_mask);
            cond1_mask = cond1_mask & trial_mask;
            cond2_mask = cond2_mask & trial_mask;
        else
            warning('Trial mask ''%s'' not found. Skipping comparison ''%s''.', ...
                    comparison.trial_mask, comp_name);
            return;
        end
    end

    %% Perform ROC Analysis
    % Check if there are enough trials in each condition
    if sum(cond1_mask) < 2 || sum(cond2_mask) < 2
        % Not enough trials to compare, so we can return empty results
        return;
    end

    % Initialize results array
    sig_vals = NaN(n_neurons, n_bins);

    % Loop through each neuron
    for i_neuron = 1:n_neurons
        % Extract data for the current neuron for both conditions
        neuron_data = squeeze(event_data.data_array(:, i_neuron, :));

        cond1_data = neuron_data(cond1_mask, :);
        cond2_data = neuron_data(cond2_mask, :);

        % Compare the two conditions for each time bin
        for i_bin = 1:n_bins
            data1 = cond1_data(:, i_bin);
            data2 = cond2_data(:, i_bin);

            % Ensure there's enough data to compare
            if sum(~isnan(data1)) > 1 && sum(~isnan(data2)) > 1
                try
                    [~, ~, ~, sig] = arrayROC(data1, data2);
                    sig_vals(i_neuron, i_bin) = sig;
                catch ME
                    % In case of error in arrayROC, leave as NaN
                end
            end
        end
    end

    %% Store Results
    analysis_results.(event_name).(comp_name).sig = sig_vals;
    analysis_results.(event_name).(comp_name).time_vector = event_data.time_vector;
    analysis_results.(event_name).(comp_name).cond_names = {cond1_name, cond2_name};

end
