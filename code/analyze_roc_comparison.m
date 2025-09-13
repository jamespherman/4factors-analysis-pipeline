%% analyze_roc_comparison.m
%
%   Performs a bin-by-bin ROC analysis to compare firing rates between two
%   specified conditions for key trial events.
%
% INPUTS:
%   core_data  - A struct containing processed neuronal data.
%   conditions - A struct with logical masks for different trial conditions.
%
% OUTPUT:
%   analysis_results - A struct containing the significance matrix for each
%                      of the three main comparisons.
%
% Author: Jules
% Date: 2025-09-09

function analysis_results = analyze_roc_comparison(core_data, conditions, is_av_session, varargin)

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Parse Optional Arguments
p = inputParser;
addParameter(p, 'comparison', struct(), @isstruct);
parse(p, varargin{:});
specific_comparison = p.Results.comparison;

%% Analysis Parameters
roc_perms = 200;
n_neurons = size(core_data.spikes.CUE_ON.rates, 1);

% Define all possible comparisons
all_comparisons = {
    {'CUE_ON', 'is_normal_dist', 'is_uniform_dist', 'Dist_at_Cue'}, ...
    {'outcomeOn', 'is_norm_common', 'is_norm_rare_high', 'RPE_at_Outcome'}, ...
    {'reward', 'is_norm_common', 'is_norm_rare_high', 'RPE_at_Reward'}, ...
    };

if is_av_session
    av_comparisons = {
        {'outcomeOn', 'is_flicker_certain', 'is_flicker_surprising', 'SPE_at_Outcome'}
        };
    all_comparisons = [all_comparisons, av_comparisons];
end

try
% Determine which comparisons to run
if ~isempty(fieldnames(specific_comparison))
    % Run only the specified comparison
    comparisons_to_run = {
        {specific_comparison.event, specific_comparison.cond1, ...
        specific_comparison.cond2, specific_comparison.name}
        };
else
    % Run all defined comparisons
    comparisons_to_run = all_comparisons;
end

catch me
    keyboard
end

%% Pre-calculate Total Workload
total_workload = 0;
total_clusters_in_session = 0;
for i_comp_dry = 1:numel(comparisons_to_run)
    comp_params_dry = comparisons_to_run{i_comp_dry};
    align_name_dry = comp_params_dry{1};
    cond1_name_dry = comp_params_dry{2};
    cond2_name_dry = comp_params_dry{3};

    time_vector_dry = core_data.spikes.(align_name_dry).time_vector;
    n_cols = numel(time_vector_dry);

    trials_cond1_dry = conditions.(cond1_name_dry);
    trials_cond2_dry = conditions.(cond2_name_dry);

    n_rows1 = sum(trials_cond1_dry);
    n_rows2 = sum(trials_cond2_dry);

    if n_rows1 >= 5 && n_rows2 >= 5
        % The workload is estimated as (rows1 + rows2) * cols per neuron
        workload_per_neuron = (n_rows1 + n_rows2) * n_cols;
        total_workload = total_workload + (workload_per_neuron * n_neurons);
        total_clusters_in_session = total_clusters_in_session + n_neurons;
    end
end

%% Initialize Tracking Variables
workload_completed = 0;
clusters_processed = 0;
master_timer = tic;
if total_clusters_in_session > 0
    fprintf('Starting analysis for %d clusters...\n', total_clusters_in_session);
end

%% Main Analysis Loop
% Iterate through each of the three comparisons
for i_comp = 1:numel(comparisons_to_run)
    comp_params = comparisons_to_run{i_comp};
    align_name = comp_params{1};
    cond1_name = comp_params{2};
    cond2_name = comp_params{3};
    result_name = comp_params{4};

    % Get data for the current alignment
    rates_all_neurons = core_data.spikes.(align_name).rates;
    time_vector = core_data.spikes.(align_name).time_vector;
    n_time_bins = numel(time_vector);

    % Get trial indices for the two conditions
    trials_cond1 = conditions.(cond1_name);
    trials_cond2 = conditions.(cond2_name);

    % Preallocate storage for significance results
    sig_results = nan(n_neurons, n_time_bins);

    % Iterate through each neuron
    for i_neuron = 1:n_neurons
        % Extract rates for the current neuron, all trials, all bins
        neuron_rates_all_trials = squeeze(rates_all_neurons(i_neuron, :, :));

        % Get the trial-by-bin matrices for each condition
        rates_cond1 = neuron_rates_all_trials(trials_cond1, :);
        rates_cond2 = neuron_rates_all_trials(trials_cond2, :);

        % Filter out silent trials for each condition independently
        is_silent_cond1 = all(rates_cond1 == 0 | isnan(rates_cond1), 2);
        rates_cond1 = rates_cond1(~is_silent_cond1, :);

        is_silent_cond2 = all(rates_cond2 == 0 | isnan(rates_cond2), 2);
        rates_cond2 = rates_cond2(~is_silent_cond2, :);

        % Skip if there are not enough trials in either condition
        if size(rates_cond1, 1) < 5 || size(rates_cond2, 1) < 5
            continue;
        end

        % Perform bin-by-bin ROC analysis
        % Note: arrayROC is efficient for this; no need to loop bins
        [~, ~, ~, sig] = arrayROC(rates_cond1, rates_cond2, roc_perms);
        sig_results(i_neuron, :) = sig;

        % Update completed workload and cluster count
        workload_per_call = (size(rates_cond1, 1) + size(rates_cond2, 1)) * size(rates_cond1, 2);
        workload_completed = workload_completed + workload_per_call;
        clusters_processed = clusters_processed + 1;

        % --- Progress Reporting ---
        elapsed_time = toc(master_timer);
        rate = workload_completed / elapsed_time;
        etc_seconds = (total_workload - workload_completed) / rate;
        etc_minutes = floor(etc_seconds / 60);
        etc_rem_seconds = rem(etc_seconds, 60);

        fprintf('Analysis completed for %d/%d clusters. Approximately %d minutes and %.0f seconds remaining.\r', ...
            clusters_processed, total_clusters_in_session, etc_minutes, etc_rem_seconds);
    end

    % Store the results in the new nested structure
    analysis_results.(align_name).(result_name).sig = sig_results;
    analysis_results.(align_name).(result_name).time_vector = time_vector;
    analysis_results.(align_name).(result_name).cond_names = ...
        {cond1_name, cond2_name};
end

if total_clusters_in_session > 0
    fprintf('\nAnalysis complete.\n');
end

end
