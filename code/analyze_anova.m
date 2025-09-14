%% analyze_anova.m
%
% Performs an N-way ANOVA on neuronal firing rates based on a provided
% analysis plan. This function is designed to be dynamically controlled by
% the settings in the `anova_plan` struct.
%
% Author: Jules
% Date: 2025-09-14
%

function session_data = analyze_anova(session_data, core_data, conditions, anova_plan)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Unpack ANOVA Parameters
event_name = anova_plan.event;
factors = anova_plan.factors;
trial_mask_name = anova_plan.trial_mask;

% Initialize results structure
session_data.analysis.anova_results = struct();

%% Prepare Data for ANOVA
% Extract data for the specified alignment event
if ~isfield(core_data, event_name)
    warning('analyze_anova:eventNotFound', ...
        'Event ''%s'' not found in core_data. Skipping ANOVA.', ...
        event_name);
    return;
end
event_data = core_data.(event_name);
[~, ~, ~] = size(event_data.binned_spikes);

% Calculate mean firing rate in the post-event window
time_vector = event_data.time_vector;
post_event_bins = time_vector >= 0;
mean_fr = mean(event_data.binned_spikes(:, :, post_event_bins), 3, 'omitnan');

% Apply the trial mask to select trials for the ANOVA
if ~isfield(conditions, trial_mask_name)
    warning('analyze_anova:maskNotFound', ...
        'Trial mask ''%s'' not found in conditions. Skipping ANOVA.', ...
        trial_mask_name);
    return;
end
trial_mask = conditions.(trial_mask_name);

% Filter the firing rates by the trial mask
filtered_fr = mean_fr(trial_mask, :);

%% Prepare Grouping Variables and Run ANOVA
% Create grouping variables for the ANOVA
group_vars = cell(1, length(factors));
for i = 1:length(factors)
    factor_name = factors{i};

    % Define the condition names based on the factor name
    high_cond_name = ['is_high_' factor_name];
    low_cond_name = ['is_low_' factor_name];

    % Check if the conditions exist
    if ~isfield(conditions, high_cond_name) || ~isfield(conditions, low_cond_name)
        warning('analyze_anova:conditionNotFound', ...
            'Condition for factor ''%s'' not found. Skipping ANOVA.', ...
            factor_name);
        return;
    end

    % Get the condition masks and filter them by the trial_mask
    high_mask = conditions.(high_cond_name)(trial_mask);
    low_mask = conditions.(low_cond_name)(trial_mask);

    % Create the grouping variable
    group_var = cell(length(high_mask), 1);
    group_var(high_mask) = {'high'};
    group_var(low_mask) = {'low'};
    group_vars{i} = group_var;
end

% Run ANOVA for each neuron
n_neurons = size(filtered_fr, 2);
p_values = cell(1, n_neurons);
term_names = {};
term_names_collected = false;

for i_neuron = 1:n_neurons
    y = filtered_fr(:, i_neuron);

    % Skip if there are not enough data points
    if sum(~isnan(y)) < length(y) * 0.5 || length(y) < 2
        continue;
    end

    [p, tbl, ~, ~] = anovan(y, group_vars, 'model', 'interaction', ...
        'varnames', factors, 'display', 'off');

    p_values{i_neuron} = p;

    if ~term_names_collected && ~isempty(tbl)
        term_names = tbl(2:(end-2), 1); % Get term names from the table
        term_names_collected = true;
    end
end

%% Store Results
if term_names_collected
    % Initialize fields in the results struct
    for i_term = 1:length(term_names)
        field_name = matlab.lang.makeValidName(term_names{i_term});
        session_data.analysis.anova_results.(field_name) = NaN(1, n_neurons);
    end

    % Populate the results struct with p-values
    for i_neuron = 1:n_neurons
        if ~isempty(p_values{i_neuron})
            for i_term = 1:length(term_names)
                field_name = matlab.lang.makeValidName(term_names{i_term});
                p_val = p_values{i_neuron}(i_term);
                session_data.analysis.anova_results.(field_name)(i_neuron) = p_val;
            end
        end
    end
end

end
