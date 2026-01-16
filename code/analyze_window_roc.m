%% analyze_window_roc.m
%
% Computes ROC AUC values in fixed time windows for each neuron,
% comparing high vs low conditions for each factor (reward, salience,
% probability, identity). This provides a single summary metric per
% neuron × factor × epoch combination, suitable for scatter plots and
% population-level comparisons.
%
% INPUTS:
%   session_data   - Standard session data structure. Must contain
%                    metadata.brain_area ('SC' or 'SNc') to select
%                    appropriate epoch windows.
%   conditions     - Output of define_task_conditions containing
%                    logical trial masks.
%   condition_defs - Output of define_task_conditions containing
%                    analysis configuration including window_roc_plan.
%   core_data      - Binned spike data from prepare_core_data. Expected
%                    structure: core_data.spikes.<event>.rates is
%                    [n_neurons × n_trials × n_time_bins] and
%                    core_data.spikes.<event>.time_vector gives bin centers.
%
% OUTPUT:
%   window_roc - Nested struct organized as window_roc.<epoch>.<factor>
%                with fields:
%                  .auc:      [n_neurons × 1] ROC AUC values (0.5 = no discrimination)
%                  .p:        [n_neurons × 1] p-values from permutation test
%                  .ci:       [n_neurons × 2] confidence intervals [lower, upper]
%                  .n_trials: [n_neurons × 2] trial counts [n_cond1, n_cond2]
%
% Author: Claude Code
% Date: 2026-01-15

function window_roc = analyze_window_roc(session_data, conditions, ...
    condition_defs, core_data)

%% Setup
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

% Get analysis configuration
wroc_plan = condition_defs.window_roc_plan;
params = wroc_plan.params;

% Determine brain area for window selection
brain_area = session_data.metadata.brain_area;
if ~ismember(brain_area, {'SC', 'SNc'})
    warning('analyze_window_roc:unknownBrainArea', ...
        'Unknown brain_area ''%s''. Defaulting to SC windows.', brain_area);
    brain_area = 'SC';
end

% Get epoch names
epoch_names = fieldnames(wroc_plan.epoch_windows);

% Get number of neurons from the first available event
first_event = wroc_plan.epoch_windows.(epoch_names{1}).event;
if ~isfield(core_data.spikes, first_event)
    error('analyze_window_roc:missingEvent', ...
        'Event ''%s'' not found in core_data.spikes.', first_event);
end
n_neurons = size(core_data.spikes.(first_event).rates, 1);

% Initialize output structure
window_roc = struct();

%% Process Each Epoch
for i_epoch = 1:length(epoch_names)
    epoch_name = epoch_names{i_epoch};
    epoch_config = wroc_plan.epoch_windows.(epoch_name);

    % Get alignment event and window for this brain area
    event_name = epoch_config.event;
    time_window = epoch_config.(brain_area);  % [start_sec, end_sec]

    % Check if event exists in core_data
    if ~isfield(core_data.spikes, event_name)
        warning('analyze_window_roc:missingEvent', ...
            'Event ''%s'' for epoch ''%s'' not found. Skipping.', ...
            event_name, epoch_name);
        continue;
    end

    event_data = core_data.spikes.(event_name);
    time_vector = event_data.time_vector;
    rates = event_data.rates;  % [n_neurons × n_trials × n_time_bins]

    % Find bins within the analysis window
    bin_mask = (time_vector >= time_window(1)) & ...
               (time_vector <= time_window(2));

    if sum(bin_mask) == 0
        warning('analyze_window_roc:noValidBins', ...
            'No bins found within window [%.3f, %.3f] for epoch ''%s''. Skipping.', ...
            time_window(1), time_window(2), epoch_name);
        continue;
    end

    % Average firing rate across bins in window for each neuron and trial
    % Result: [n_neurons × n_trials]
    window_rates = mean(rates(:, :, bin_mask), 3, 'omitnan');

    %% Process Each Factor
    for i_factor = 1:length(wroc_plan.factors)
        factor_config = wroc_plan.factors(i_factor);
        factor_name = factor_config.name;

        % Initialize results for this factor
        auc_vals = NaN(n_neurons, 1);
        p_vals = NaN(n_neurons, 1);
        ci_vals = NaN(n_neurons, 2);
        n_trials_vals = NaN(n_neurons, 2);

        % Get condition masks
        if ~isfield(conditions, factor_config.cond1) || ...
           ~isfield(conditions, factor_config.cond2)
            warning('analyze_window_roc:missingCondition', ...
                'Condition masks for factor ''%s'' not found. Skipping.', ...
                factor_name);
            continue;
        end

        cond1_mask = conditions.(factor_config.cond1);
        cond2_mask = conditions.(factor_config.cond2);

        % Apply trial subset masks (AND all masks together)
        trial_masks = factor_config.trial_mask;
        if ~isempty(trial_masks)
            for i_mask = 1:length(trial_masks)
                mask_name = trial_masks{i_mask};
                if isfield(conditions, mask_name)
                    cond1_mask = cond1_mask & conditions.(mask_name);
                    cond2_mask = cond2_mask & conditions.(mask_name);
                else
                    warning('analyze_window_roc:missingTrialMask', ...
                        'Trial mask ''%s'' not found for factor ''%s''.', ...
                        mask_name, factor_name);
                end
            end
        end

        %% Process Each Neuron
        for i_neuron = 1:n_neurons
            % Get window-averaged firing rates for this neuron
            neuron_rates = window_rates(i_neuron, :)';  % [n_trials × 1]

            % Extract data for each condition
            cond1_data = neuron_rates(cond1_mask);
            cond2_data = neuron_rates(cond2_mask);

            % Remove NaN trials (trials with missing event times)
            cond1_data = cond1_data(~isnan(cond1_data));
            cond2_data = cond2_data(~isnan(cond2_data));

            % Store trial counts
            n1 = length(cond1_data);
            n2 = length(cond2_data);
            n_trials_vals(i_neuron, :) = [n1, n2];

            % Check minimum trial requirement
            if n1 < params.min_trials || n2 < params.min_trials
                % Leave as NaN
                continue;
            end

            % Compute ROC AUC using arrayROC
            % Note: arrayROC expects [nTrials × nBins], so we transpose
            % to get [n1 × 1] and [n2 × 1] which works as single-bin data
            [auc, ci, ~, ~] = arrayROC(cond1_data, cond2_data, ...
                params.n_bootstrap, params.alpha, 'notpfp');

            auc_vals(i_neuron) = auc;
            ci_vals(i_neuron, :) = ci';  % ci is [2 × 1], transpose to [1 × 2]

            % Compute p-value using permutation test
            p_vals(i_neuron) = permutation_p_value(cond1_data, cond2_data, ...
                auc, params.n_bootstrap);
        end

        % Store results for this factor
        window_roc.(epoch_name).(factor_name).auc = auc_vals;
        window_roc.(epoch_name).(factor_name).p = p_vals;
        window_roc.(epoch_name).(factor_name).ci = ci_vals;
        window_roc.(epoch_name).(factor_name).n_trials = n_trials_vals;
    end
end

end

%% ========================================================================
%  PERMUTATION P-VALUE COMPUTATION
%  ========================================================================
function p = permutation_p_value(x1, x2, observed_auc, n_perm)
% Compute two-sided p-value using permutation test
%
% Inputs:
%   x1, x2       - Data vectors for each condition
%   observed_auc - The observed AUC value
%   n_perm       - Number of permutations
%
% Output:
%   p - Two-sided p-value

% Combine data
pooled = [x1; x2];
n1 = length(x1);
n_total = length(pooled);

% Generate permuted AUCs
perm_aucs = zeros(n_perm, 1);
for i = 1:n_perm
    % Shuffle the pooled data
    perm_idx = randperm(n_total);
    perm1 = pooled(perm_idx(1:n1));
    perm2 = pooled(perm_idx(n1+1:end));

    % Compute AUC for permuted data
    perm_aucs(i) = compute_auc_fast(perm1, perm2);
end

% Compute two-sided p-value
% Count how many permuted AUCs are as extreme as or more extreme than observed
% "Extreme" means far from 0.5 in either direction
observed_deviation = abs(observed_auc - 0.5);
perm_deviations = abs(perm_aucs - 0.5);

p = (sum(perm_deviations >= observed_deviation) + 1) / (n_perm + 1);

end

%% ========================================================================
%  FAST AUC COMPUTATION (WITHOUT CI)
%  ========================================================================
function auc = compute_auc_fast(x1, x2)
% Compute ROC AUC quickly without bootstrap CI
% Uses the Mann-Whitney U statistic relationship to AUC
%
% AUC = P(X2 > X1) = U / (n1 * n2)
% where U is the Mann-Whitney U statistic

n1 = length(x1);
n2 = length(x2);

if n1 < 2 || n2 < 2
    auc = NaN;
    return;
end

% Compute using ranks (more efficient for larger samples)
% Pool data and rank
pooled = [x1(:); x2(:)];
[~, sort_idx] = sort(pooled);
ranks = zeros(size(pooled));
ranks(sort_idx) = 1:length(pooled);

% Sum of ranks for condition 2
rank_sum_2 = sum(ranks(n1+1:end));

% Mann-Whitney U for condition 2 (tests if x2 > x1)
U2 = rank_sum_2 - n2*(n2+1)/2;

% AUC = U / (n1 * n2)
auc = U2 / (n1 * n2);

end
