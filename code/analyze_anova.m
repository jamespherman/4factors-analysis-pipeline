%% analyze_anova.m
%
% Performs a time-resolved 4-way ANOVA on neuronal firing rates for the
% 4-factors task. This function iterates through each time bin of the
% aligned neural data, performing a separate ANOVA at each point.
%
% This function uses a fixed 4-way interaction model with factors
% determined by the `anova_plan`. It retrieves pre-computed categorical
% factors from the `conditions.factors` struct.
%
% Inputs:
%   session_data - The main data struct for the session.
%   core_data    - The core data structure containing aligned neural data.
%   conditions   - A struct with logical masks for trial conditions and
%                  a 'factors' sub-struct with categorical variables.
%   anova_plan   - A struct defining ANOVA parameters, including the
%                  alignment event, trial mask, and factors to test.
%                  The 'trial_mask' can be a string or a cell array of
%                  strings to combine multiple conditions.
%
% Output:
%   session_data - Updated struct with time-resolved ANOVA results.
%                  For each factor, results are stored as a
%                  [n_neurons x n_time_bins] matrix.
%
% Author: Jules
% Date: 2025-09-21 (Refactored for time-resolved analysis)
%

function session_data = analyze_anova(session_data, core_data, ...
    conditions, anova_plan)

    %% Setup Paths
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'utils'));

    %% Unpack ANOVA Parameters
    event_name = anova_plan.event;
    analysis_name = anova_plan.name; % e.g., 'anova_imagetrials'

    %% Initialize Results Structure
    if ~isfield(session_data, 'analysis')
        session_data.analysis = struct();
    end
    if ~isfield(session_data.analysis, 'anova_results')
        session_data.analysis.anova_results = struct();
    end
    session_data.analysis.anova_results.(analysis_name).(event_name) = ...
        struct();

    %% Prepare Data for ANOVA
    if ~isfield(core_data.spikes, event_name)
        warning('analyze_anova:eventNotFound', ...
            'Event ''%s'' not found in core_data.spikes. Skipping.', ...
            event_name);
        return;
    end
    event_data = core_data.spikes.(event_name);
    time_vector = event_data.time_vector;
    n_bins = length(time_vector);

    % Get the trial mask to select trials for the ANOVA
    trial_mask_spec = anova_plan.trial_mask;
    if ischar(trial_mask_spec)
        trial_mask_spec = {trial_mask_spec};
    end

    if ~iscell(trial_mask_spec)
        error('analyze_anova:invalidTrialMask', ['trial_mask must ' ...
            'be a string or a cell array of strings.']);
    end

    if isempty(trial_mask_spec)
        warning('analyze_anova:emptyTrialMask', ['Trial mask spec is ' ...
            '            empty. Skipping.']);
        return;
    end

    first_mask_name = trial_mask_spec{1};
    if ~isfield(conditions, first_mask_name)
        warning('analyze_anova:maskNotFound', ['Trial mask ''%s'' ' ...
            '            not found. Skipping.'], ...
            first_mask_name);
        return;
    end
    trial_mask = conditions.(first_mask_name);

    for i = 2:length(trial_mask_spec)
        mask_name = trial_mask_spec{i};
        if ~isfield(conditions, mask_name)
            warning('analyze_anova:maskNotFound', ['Trial mask ''%s'' ' ...
                '            not found. Skipping.'], ...
                mask_name);
            return;
        end
        trial_mask = trial_mask & conditions.(mask_name);
    end

    % Filter firing rates by the trial mask
    % This results in a [n_neurons x n_valid_trials x n_time_bins] matrix
    fr_data = event_data.rates(:, trial_mask, :);
    n_neurons = size(fr_data, 1);

    %% Prepare Factors for ANOVA
    factor_names_lower = anova_plan.factors;
    num_factors = length(factor_names_lower);
    group_vars = cell(1, num_factors);
    factor_names_capitalized = cell(1, num_factors);

    for i = 1:num_factors
        factor_name = factor_names_lower{i};
        if ~isfield(conditions.factors, factor_name)
            warning('analyze_anova:factorNotFound', ...
                'Factor ''%s'' not found. Skipping.', factor_name);
            return;
        end
        factor_data = conditions.factors.(factor_name);
        group_vars{i} = factor_data(trial_mask);
        capitalized_name = [upper(factor_name(1)), factor_name(2:end)];
        factor_names_capitalized{i} = capitalized_name;
    end

    %% Create a master mask for valid trials
    % Ensures only trials with complete data across all factors are used.
    valid_trials_mask = true(sum(trial_mask), 1);
    for i = 1:num_factors
        valid_trials_mask = valid_trials_mask & ~cellfun('isempty', ...
            group_vars{i});
    end

    % Filter firing rates and group variables by the valid trials mask
    fr_filt = fr_data(:, valid_trials_mask, :);
    for i = 1:num_factors
        group_vars{i} = group_vars{i}(valid_trials_mask);
    end

    if size(fr_filt, 2) < num_factors * 2
        warning('analyze_anova:notEnoughData', ['Not enough valid ' ...
            'trials. Skipping.']);
        return;
    end

    %% Initialize result matrices
    term_names = {};
    p_val_fields = {};
    f_val_fields = {};

    % Run once to get term names
    y_template = squeeze(fr_filt(1, :, 1))';
    [~, tbl, ~] = anovan(y_template, group_vars, 'model', ...
        'interaction', 'varnames', factor_names_capitalized, ...
        'display', 'off');

    if ~isempty(tbl)
        term_names = tbl(2:(end-2), 1);
        for i_term = 1:length(term_names)
            base_name = matlab.lang.makeValidName(term_names{i_term});
            p_val_fields{i_term} = ['p_', base_name];
            f_val_fields{i_term} = ['f_', base_name];

            % Initialize result matrices with NaNs
            session_data.analysis.anova_results.(analysis_name).(...
                event_name).(p_val_fields{i_term}) = NaN(n_neurons, ...
                n_bins);
            session_data.analysis.anova_results.(analysis_name).(...
                event_name).(f_val_fields{i_term}) = NaN(n_neurons, ...
                n_bins);
        end
    end

    if isempty(term_names)
        warning('analyze_anova:noTerms', ['Could not determine ANOVA ' ...
            '            terms. Skipping.']);
        return;
    end

    %% Time-Resolved ANOVA Loop
    for i_bin = 1:n_bins
        for i_neuron = 1:n_neurons
            y = squeeze(fr_filt(i_neuron, :, i_bin))';

            if sum(~isnan(y)) < length(y) * 0.5 || isempty(y) || ...
                    all(y == y(1))
                continue; % Skip if data is insufficient or constant
            end

            [p, tbl, ~, ~] = anovan(y, group_vars, 'model', ...
                'interaction', 'varnames', factor_names_capitalized, ...
                'display', 'off');

            % Store p-values and F-statistics
            for i_term = 1:length(term_names)
                p_val = p(i_term);
                f_val = tbl{i_term + 1, 6}; % F-statistic is in column 6

                session_data.analysis.anova_results.(analysis_name).(...
                    event_name).(p_val_fields{i_term})(i_neuron, i_bin)...
                    = p_val;
                session_data.analysis.anova_results.(analysis_name).(...
                    event_name).(f_val_fields{i_term})(i_neuron, i_bin)...
                    = f_val;
            end
        end
    end

    % Store the time vector for plotting
    session_data.analysis.anova_results.(analysis_name).(...
        event_name).time_vector = time_vector;
end