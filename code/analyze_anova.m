%% analyze_anova.m
%
% Performs an N-way ANOVA on neuronal firing rates. This function is
% dynamically driven by the `anova_plan` struct, which specifies the
% alignment event, factors, and trials to include.
%
% Inputs:
%   session_data - The main data struct for the session, which will be
%                  updated with the analysis results.
%   core_data    - The core data structure containing aligned neural data.
%   conditions   - A struct with logical masks for trial conditions.
%   anova_plan   - A struct defining the ANOVA to be run. It must contain:
%                  .event: The name of the alignment event (e.g., 'CUE_ON').
%                  .trial_mask: The name of the condition mask to select
%                               trials for the analysis.
%                  .factors: A nested cell array defining the factors. Each
%                            inner cell array should contain the full names
%                            of the condition masks that represent the
%                            levels of a single factor.
%                            Example: {{'cond_A1', 'cond_A2'}, {'cond_B1', 'cond_B2'}}
%
% Output:
%   session_data - The input session_data struct, with the ANOVA results
%                  added under the `analysis.anova_results.(event_name)` field.
%
% Author: Jules
% Date: 2025-09-19
%

function session_data = analyze_anova(session_data, core_data, ...
    conditions, anova_plan)

    %% Setup Paths
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'utils'));

    %% Unpack ANOVA Parameters
    event_name = anova_plan.event;
    factors = anova_plan.factors;
    trial_mask_name = anova_plan.trial_mask;

    %% Initialize Results Structure
    if ~isfield(session_data, 'analysis')
        session_data.analysis = struct();
    end
    session_data.analysis.anova_results.(event_name) = struct();

    %% Prepare Data for ANOVA
    if ~isfield(core_data.spikes, event_name)
        warning('analyze_anova:eventNotFound', ...
            'Event ''%s'' not found in core_data. Skipping ANOVA.', ...
            event_name);
        return;
    end
    event_data = core_data.spikes.(event_name);
    time_vector = event_data.time_vector;

    % Calculate mean firing rate in the post-event window (t >= 0)
    post_event_bins = time_vector >= 0;
    mean_fr = mean(event_data.rates(:, :, post_event_bins), 3, 'omitnan');

    % Get the trial mask to select trials for the ANOVA
    if ~isfield(conditions, trial_mask_name)
        warning('analyze_anova:maskNotFound', ...
            'Trial mask ''%s'' not found in conditions. Skipping ANOVA.', ...
            trial_mask_name);
        return;
    end
    trial_mask = conditions.(trial_mask_name);
    filtered_fr = mean_fr(:, trial_mask);

    %% Prepare Grouping Variables and Run ANOVA
    group_vars = cell(1, length(factors));
    factor_names = cell(1, length(factors));

    for i = 1:length(factors)
        factor_levels = factors{i};
        % Dynamically derive a factor name (e.g., 'cond_A1' -> 'condA')
        clean_name = strrep(factor_levels{1}, '_', '');
        factor_names{i} = matlab.lang.makeValidName(clean_name(1:end-1));

        group_var = cell(sum(trial_mask), 1);
        for j = 1:length(factor_levels)
            level_name = factor_levels{j};
            if ~isfield(conditions, level_name)
                warning('analyze_anova:levelNotFound', ...
                    'Condition ''%s'' not found. Skipping ANOVA.', level_name);
                return;
            end
            level_mask = conditions.(level_name)(trial_mask);
            group_var(level_mask) = {level_name};
        end
        group_vars{i} = group_var;
    end

    n_neurons = size(filtered_fr, 1);
    p_values = cell(1, n_neurons);
    term_names_collected = false;
    tbl_term_names = {};

    for i_neuron = 1:n_neurons
        y = filtered_fr(i_neuron, :)';
        if sum(~isnan(y)) < length(y) * 0.5 || isempty(y)
            continue;
        end

        [p, tbl, ~, ~] = anovan(y, group_vars, 'model', 'interaction', ...
            'varnames', factor_names, 'display', 'off');
        p_values{i_neuron} = p;

        if ~term_names_collected && ~isempty(tbl)
            tbl_term_names = tbl(2:(end-2), 1);
            term_names_collected = true;
        end
    end

    %% Store Results
    if term_names_collected
        for i_term = 1:length(tbl_term_names)
            field_name = matlab.lang.makeValidName(tbl_term_names{i_term});
            session_data.analysis.anova_results.(event_name).(field_name) = NaN(1, n_neurons);
        end

        for i_neuron = 1:n_neurons
            if ~isempty(p_values{i_neuron})
                for i_term = 1:length(tbl_term_names)
                    field_name = matlab.lang.makeValidName(tbl_term_names{i_term});
                    p_val = p_values{i_neuron}(i_term);
                    session_data.analysis.anova_results.(event_name).(field_name)(i_neuron) = p_val;
                end
            end
        end
    end
end
