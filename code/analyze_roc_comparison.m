%% analyze_roc_comparison.m
%
% This function performs a plan-driven ROC analysis to compare neural
% firing rates between two specified conditions.
%
% Inputs:
%   core_data       - The core data structure containing aligned neural data.
%   conditions      - A struct containing logical masks for various trial
%                     conditions.
%   roc_plan_item   - A struct from the analysis plan detailing the
%                     comparison, with fields: event, cond1, cond2, name,
%                     and an optional trial_mask.
%
% Output:
%   analysis_results - A struct containing the ROC analysis results,
%                      including p-values (sig), a time vector, and the
%                      names of the conditions that were compared.
%
% Author: Jules
% Date: 2025-09-23
%

function analysis_results = analyze_roc_comparison(core_data, ...
    conditions, roc_plan_item)
    %% Setup Paths
    % Add the 'utils' directory to the path so that helper functions can be
    % found.
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'utils'));

    %% Initialize Output
    analysis_results = struct();

    %% Extract Comparison Parameters
    event_name = roc_plan_item.event;
    cond1_name = roc_plan_item.cond1;
    cond2_name = roc_plan_item.cond2;
    comp_name = roc_plan_item.name;

    % Check if the alignment event exists in the core data
    if ~isfield(core_data.spikes, event_name)
        warning('Event ''%s'' not found in core_data.spikes. Skipping ''%s''.', ...
                event_name, comp_name);
        return;
    end

    event_data = core_data.spikes.(event_name);
    [n_neurons, ~, n_bins] = size(event_data.rates);

    %% Get Trial Masks
    if ~isfield(conditions, cond1_name) || ~isfield(conditions, cond2_name)
        warning('Conds ''%s'' or ''%s'' not found. Skipping ''%s''.', ...
                cond1_name, cond2_name, comp_name);
        return;
    end

    cond1_mask = conditions.(cond1_name);
    cond2_mask = conditions.(cond2_name);

    % Apply an overall trial mask if specified
    if isfield(roc_plan_item, 'trial_mask') && ~isempty(...
            roc_plan_item.trial_mask)
        if isfield(conditions, roc_plan_item.trial_mask)
            trial_mask = conditions.(roc_plan_item.trial_mask);
            cond1_mask = cond1_mask & trial_mask;
            cond2_mask = cond2_mask & trial_mask;
        else
            warning('Trial mask ''%s'' not found. Skipping ''%s''.', ...
                    roc_plan_item.trial_mask, comp_name);
            return;
        end
    end

    %% Perform ROC Analysis
    if sum(cond1_mask) < 2 || sum(cond2_mask) < 2
        return;
    end

    sig_vals = NaN(n_neurons, n_bins);

    for i_neuron = 1:n_neurons
        neuron_data = squeeze(event_data.rates(i_neuron, :, :));
        cond1_data = neuron_data(cond1_mask, :);
        cond2_data = neuron_data(cond2_mask, :);

        for i_bin = 1:n_bins
            data1 = cond1_data(:, i_bin);
            data2 = cond2_data(:, i_bin);

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
    analysis_results.sig = sig_vals;
    analysis_results.time_vector = event_data.time_vector;
    analysis_results.cond_names = {cond1_name, cond2_name};

end
