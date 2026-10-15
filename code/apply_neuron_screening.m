function [selected_neurons, screening_info, session_data] = apply_neuron_screening(session_data, core_data, condition_defs, project_root)
% APPLY_NEURON_SCREENING Apply refined neuron screening criteria.
%
%   [selected_neurons, screening_info] = APPLY_NEURON_SCREENING(session_data,
%       core_data, condition_defs)
%   [selected_neurons, screening_info, session_data] = APPLY_NEURON_SCREENING(session_data,
%       core_data, condition_defs, project_root)
%
%   This function implements the three-part screening criteria:
%     1. Mean firing rate threshold (brain-area specific, computed over active period)
%     2. Task modulation (significant deviation from baseline at ANY location/event)
%     3. Sparsity threshold (proportion of empty 100ms bins must be < 0.7)
%
%   A neuron is included if ALL THREE conditions are met.
%
%   Active Period: Metrics are computed over each neuron's active period
%   (first spike to last spike), not the full recording duration. This avoids
%   penalizing neurons that were lost partway through or only became isolated later.
%
%   When use_strict_screening is false, all neurons are marked as included,
%   but screening_info is still computed for all neurons.
%
%   For SC sessions, this function also calls screen_sc_neurons to obtain the
%   detailed neuron_class categorization ('classic', 'interneuron', 'excluded').
%
%   CACHING: This function caches modulation_info in session_data.cached_modulation_info
%   to avoid expensive recomputation on subsequent runs. The cache is validated
%   by checking that cluster_ids match.
%
%   INPUTS:
%       session_data   - Standard session data structure (must include
%                        spikes.times, spikes.clusters, trialInfo, metadata.brain_area)
%       core_data      - Binned spike data from prepare_core_data
%       condition_defs - Analysis configuration from define_task_conditions
%       project_root   - (Optional) Project root path, needed for SC classification
%
%   OUTPUTS:
%       selected_neurons - Logical vector (nNeurons x 1) where true means included
%       screening_info   - Struct array (one per neuron) with fields:
%       session_data     - Updated session_data with cached_modulation_info (for saving)
%           .cluster_id              - Cluster ID
%           .included                - true/false
%           .neuron_class            - 'classic', 'interneuron', 'excluded', or 'modulated' (SNc)
%           .mean_firing_rate        - Mean FR over active period (sp/s)
%           .proportion_empty_bins   - Proportion of 100ms bins with zero spikes (0-1)
%           .passes_fr_threshold     - true/false
%           .passes_sparsity_threshold - true/false
%           .is_modulated            - true/false
%           .modulated_at            - Cell array describing where modulation detected
%           .exclusion_reason        - String describing which criteria failed
%           .active_start            - First spike time (seconds)
%           .active_end              - Last spike time (seconds)
%           .active_duration         - Duration of active period (seconds)
%           .fr_threshold            - Threshold used for this neuron
%           .sparsity_threshold      - Sparsity threshold used
%           .recording_duration      - Total recording duration (for reference)
%
%   Author: Claude Code
%   Date: 2026-01-18

fprintf('apply_neuron_screening: Starting refined screening...\n');

%% Handle optional project_root parameter
if nargin < 4
    project_root = '';
end

%% Extract configuration
neuron_inclusion = condition_defs.neuron_inclusion;
use_strict_screening = neuron_inclusion.use_strict_screening;

% Get brain area to determine FR threshold
brain_area = session_data.metadata.brain_area;
if strcmp(brain_area, 'SC')
    fr_threshold = neuron_inclusion.min_firing_rate_sc;
elseif strcmp(brain_area, 'SNc')
    fr_threshold = neuron_inclusion.min_firing_rate_snc;
else
    warning('apply_neuron_screening: Unknown brain_area "%s". Using SC threshold.', brain_area);
    fr_threshold = neuron_inclusion.min_firing_rate_sc;
end

% Get sparsity threshold
max_proportion_empty = neuron_inclusion.max_proportion_empty_bins;

%% Get neuron information
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nNeurons = numel(cluster_ids);
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

fprintf('apply_neuron_screening: Brain area=%s, FR threshold=%.1f sp/s, max empty bins=%.0f%%, %d neurons\n', ...
    brain_area, fr_threshold, max_proportion_empty * 100, nNeurons);

%% Detect stable trial periods for each neuron (Phase 2 improvement)
% Get stability detection parameters (use defaults if not specified)
if isfield(neuron_inclusion, 'stability_detection')
    stability_params = neuron_inclusion.stability_detection;
else
    stability_params = struct(...
        'window_size', 20, ...
        'step_size', 5, ...
        'cv_threshold', 0.5, ...
        'min_fr_threshold', 0.1, ...
        'min_stable_trials', 50);
end
if ~isfield(stability_params, 'trial_scope')
    stability_params.trial_scope = 'all';
end
if ~isfield(stability_params, 'min_coverage')
    stability_params.min_coverage = 0;
end

% Trials over which stability is assessed (session trial indices)
scope_trials = get_stability_scope_trials(session_data, ...
    stability_params.trial_scope);

% Check for cached stability_info first. The cache is only valid if it
% was computed for the same neurons and the same trial scope.
use_cached_stability = false;
if isfield(session_data, 'stability_info') && ...
        ~isempty(session_data.stability_info)
    cached_stability = session_data.stability_info;
    if length(cached_stability) == nNeurons && ...
            isfield(cached_stability, 'trial_scope') && ...
            strcmp(cached_stability(1).trial_scope, ...
            stability_params.trial_scope)
        cached_ids = [cached_stability.cluster_id];
        if isequal(cached_ids(:), cluster_ids(:))
            use_cached_stability = true;
            fprintf('apply_neuron_screening: Using cached stability_info.\n');
            stability_info = cached_stability;
        end
    end
end

if ~use_cached_stability
    fprintf(['apply_neuron_screening: Detecting stable trial ' ...
        'periods over %d trials (scope: %s)...\n'], ...
        numel(scope_trials), stability_params.trial_scope);

    % Extract timing of the in-scope trials
    trial_start_times = session_data.eventTimes.trialBegin(scope_trials);
    trial_end_times = session_data.eventTimes.trialEnd(scope_trials);

    % Call the stability detection function
    stability_info = detect_stable_trial_period(...
        all_spike_times, ...
        all_spike_clusters, ...
        trial_start_times(:), ...
        trial_end_times(:), ...
        cluster_ids, ...
        stability_params);

    % Map stable periods from in-scope trial positions back to session
    % trial indices (downstream code indexes session trials), and record
    % the fraction of in-scope trials each stable period covers.
    n_scope = numel(scope_trials);
    for i = 1:nNeurons
        stability_info(i).trial_scope = stability_params.trial_scope;
        if stability_info(i).is_stable
            stability_info(i).coverage = ...
                stability_info(i).n_stable_trials / n_scope;
        else
            stability_info(i).coverage = 0;
        end
        stability_info(i).stable_trial_start = ...
            scope_trials(stability_info(i).stable_trial_start);
        stability_info(i).stable_trial_end = ...
            scope_trials(stability_info(i).stable_trial_end);
    end

    % Cache for future runs. Modulation was computed using the previous
    % stable periods, so its cache is no longer valid.
    session_data.stability_info = stability_info;
    if isfield(session_data, 'cached_modulation_info')
        session_data = rmfield(session_data, 'cached_modulation_info');
    end
end

%% Compute total recording duration from trial times (for reference)
recording_start = session_data.eventTimes.trialBegin(1);
recording_end = session_data.eventTimes.trialEnd(end);
recording_duration = recording_end - recording_start;

fprintf('apply_neuron_screening: Total recording duration = %.1f s (from trial times)\n', recording_duration);

%% Initialize outputs
selected_neurons = false(nNeurons, 1);
screening_info = struct();

for i = 1:nNeurons
    screening_info(i).cluster_id = cluster_ids(i);
    screening_info(i).included = false;
    screening_info(i).neuron_class = 'excluded';  % Default, will be updated
    screening_info(i).mean_firing_rate = 0;
    screening_info(i).proportion_empty_bins = 1;  % Default to all empty
    screening_info(i).passes_fr_threshold = false;
    screening_info(i).passes_sparsity_threshold = false;
    screening_info(i).is_modulated = false;
    screening_info(i).modulated_at = {};
    screening_info(i).modulation_details = struct();
    screening_info(i).exclusion_reason = '';
    screening_info(i).active_start = NaN;
    screening_info(i).active_end = NaN;
    screening_info(i).active_duration = 0;
    screening_info(i).fr_threshold = fr_threshold;
    screening_info(i).sparsity_threshold = max_proportion_empty;
    screening_info(i).recording_duration = recording_duration;
    screening_info(i).used_stable_period = false;  % Phase 2: track if stable period was used
    screening_info(i).stable_trial_range = [];     % Phase 2: [start_trial, end_trial]
    % v2 data-driven classification (SC only)
    screening_info(i).neuron_class_v2 = '';
end

%% Criterion 1 & 3: Mean Firing Rate and Sparsity (computed over STABLE or ACTIVE PERIOD)
% Phase 2 improvement: Use stable period if detected, otherwise fall back to active period
for i = 1:nNeurons
    cluster_id = cluster_ids(i);
    spike_times = all_spike_times(all_spike_clusters == cluster_id);
    n_spikes = length(spike_times);

    if n_spikes < 2
        % Not enough spikes to define an active period
        screening_info(i).mean_firing_rate = 0;
        screening_info(i).proportion_empty_bins = 1;
        screening_info(i).passes_fr_threshold = false;
        screening_info(i).passes_sparsity_threshold = false;
        screening_info(i).active_start = NaN;
        screening_info(i).active_end = NaN;
        screening_info(i).active_duration = 0;
        screening_info(i).used_stable_period = false;
        screening_info(i).stable_trial_range = [];
        continue;
    end

    % Check if stable period is available (Phase 2 improvement)
    if exist('stability_info', 'var') && stability_info(i).is_stable
        % Use stable period (trial-based timing)
        stable_trial_start = stability_info(i).stable_trial_start;
        stable_trial_end = stability_info(i).stable_trial_end;
        active_start = session_data.eventTimes.trialBegin(stable_trial_start);
        active_end = session_data.eventTimes.trialEnd(stable_trial_end);
        screening_info(i).used_stable_period = true;
        screening_info(i).stable_trial_range = [stable_trial_start, stable_trial_end];

        % Filter spikes to stable period
        spike_times_filtered = spike_times(spike_times >= active_start & spike_times < active_end);
        n_spikes_filtered = length(spike_times_filtered);
    else
        % Fall back to first-to-last spike (original behavior)
        active_start = spike_times(1);
        active_end = spike_times(end);
        screening_info(i).used_stable_period = false;
        screening_info(i).stable_trial_range = [];
        spike_times_filtered = spike_times;
        n_spikes_filtered = n_spikes;
    end

    active_duration = active_end - active_start;

    % Store active period info
    screening_info(i).active_start = active_start;
    screening_info(i).active_end = active_end;
    screening_info(i).active_duration = active_duration;

    % Handle edge case: very short active period
    if active_duration < 1  % Less than 1 second
        screening_info(i).mean_firing_rate = n_spikes_filtered;  % Treat as if 1 second
        screening_info(i).proportion_empty_bins = 0;    % Can't compute meaningfully
        screening_info(i).passes_fr_threshold = n_spikes_filtered >= fr_threshold;
        screening_info(i).passes_sparsity_threshold = true;  % Benefit of doubt
        continue;
    end

    % Compute mean firing rate over active/stable period only
    mean_fr = n_spikes_filtered / active_duration;
    screening_info(i).mean_firing_rate = mean_fr;
    screening_info(i).passes_fr_threshold = mean_fr >= fr_threshold;

    % Compute sparsity over active/stable period only (100ms bins)
    edges = active_start : 0.1 : active_end;
    if length(edges) < 2
        % Not enough bins
        screening_info(i).proportion_empty_bins = 0;
        screening_info(i).passes_sparsity_threshold = true;
    else
        counts = histcounts(spike_times_filtered, edges);
        proportion_empty = sum(counts == 0) / numel(counts);
        screening_info(i).proportion_empty_bins = proportion_empty;
        screening_info(i).passes_sparsity_threshold = proportion_empty < max_proportion_empty;
    end
end

n_pass_fr = sum([screening_info.passes_fr_threshold]);
n_pass_sparsity = sum([screening_info.passes_sparsity_threshold]);
n_used_stable = sum([screening_info.used_stable_period]);
fprintf('apply_neuron_screening: %d/%d neurons (%.1f%%) pass FR threshold (mean FR >= %.1f sp/s).\n', ...
    n_pass_fr, nNeurons, 100 * n_pass_fr / nNeurons, fr_threshold);
fprintf('apply_neuron_screening: %d/%d neurons (%.1f%%) pass sparsity threshold (<%.0f%% empty bins).\n', ...
    n_pass_sparsity, nNeurons, 100 * n_pass_sparsity / nNeurons, max_proportion_empty * 100);
fprintf('apply_neuron_screening: %d/%d neurons (%.1f%%) used stable period (vs active period).\n', ...
    n_used_stable, nNeurons, 100 * n_used_stable / nNeurons);

%% Criterion 2: Task Modulation
% Check for cached modulation_info to avoid expensive recomputation
% The cache is valid if:
%   1. session_data.cached_modulation_info exists
%   2. It has the same number of neurons as current cluster_ids
%   3. The cluster_ids match (same neurons, same order)
use_cached = false;
if isfield(session_data, 'cached_modulation_info') && ~isempty(session_data.cached_modulation_info)
    cached = session_data.cached_modulation_info;
    if length(cached) == nNeurons
        % Check if cluster_ids match
        cached_ids = [cached.cluster_id];
        if isequal(cached_ids(:), cluster_ids(:))
            use_cached = true;
            fprintf('apply_neuron_screening: Using cached modulation_info (skipping recomputation).\n');
        else
            fprintf('apply_neuron_screening: Cached modulation_info cluster_ids mismatch, recomputing...\n');
        end
    else
        fprintf('apply_neuron_screening: Cached modulation_info size mismatch (%d vs %d), recomputing...\n', ...
            length(cached), nNeurons);
    end
end

if use_cached
    modulation_info = session_data.cached_modulation_info;
else
    % Pass stability_info to use neuron present trials when computing modulation
    modulation_info = compute_task_modulation(session_data, core_data, condition_defs, stability_info);
    % Store in session_data for future runs (will be saved by run_4factors_analysis)
    session_data.cached_modulation_info = modulation_info;
end

for i = 1:nNeurons
    screening_info(i).is_modulated = modulation_info(i).is_modulated;
    screening_info(i).modulated_at = modulation_info(i).modulated_at;
    screening_info(i).modulation_details = modulation_info(i).modulation_details;

    % Copy bin_significance if available
    if isfield(modulation_info(i), 'bin_significance')
        screening_info(i).bin_significance = modulation_info(i).bin_significance;
    end
end

n_modulated = sum([modulation_info.is_modulated]);
fprintf('apply_neuron_screening: %d/%d neurons (%.1f%%) show task modulation.\n', ...
    n_modulated, nNeurons, 100 * n_modulated / nNeurons);

%% Combined Screening Logic (Three-Part Test)
if use_strict_screening
    % Apply all three criteria: FR threshold AND task modulation AND sparsity
    for i = 1:nNeurons
        pass_fr = screening_info(i).passes_fr_threshold;
        is_mod = screening_info(i).is_modulated;
        pass_sparse = screening_info(i).passes_sparsity_threshold;

        if pass_fr && is_mod && pass_sparse
            % All three criteria met
            screening_info(i).included = true;
            screening_info(i).exclusion_reason = '';
        else
            % Build exclusion reason dynamically based on which criteria failed
            screening_info(i).included = false;
            failed_criteria = {};

            if ~pass_fr
                failed_criteria{end+1} = 'Low FR';
            end
            if ~is_mod
                failed_criteria{end+1} = 'No mod.';
            end
            if ~pass_sparse
                failed_criteria{end+1} = 'Sparse';
            end

            % Join with ' & ' separator
            screening_info(i).exclusion_reason = strjoin(failed_criteria, ' & ');
        end

        selected_neurons(i) = screening_info(i).included;
    end
else
    % Include all neurons when use_strict_screening is false
    fprintf('apply_neuron_screening: use_strict_screening=false, including all neurons.\n');
    for i = 1:nNeurons
        screening_info(i).included = true;
        screening_info(i).exclusion_reason = '';
    end
    selected_neurons = true(nNeurons, 1);
end

%% Summary
n_included = sum(selected_neurons);
fprintf('apply_neuron_screening: Final selection: %d/%d neurons (%.1f%%) included.\n', ...
    n_included, nNeurons, 100 * n_included / nNeurons);

% Count exclusion reasons
if use_strict_screening
    exclusion_reasons = {screening_info.exclusion_reason};
    unique_reasons = unique(exclusion_reasons(~cellfun(@isempty, exclusion_reasons)));

    fprintf('  Exclusion breakdown:\n');
    for i = 1:length(unique_reasons)
        reason = unique_reasons{i};
        count = sum(strcmp(exclusion_reasons, reason));
        fprintf('    %s: %d\n', reason, count);
    end
end

%% Brain Area-Specific Classification (neuron_class)
% For SC sessions, call screen_sc_neurons to get detailed classification:
%   'classic'     - Classic SC neuron (contralateral, not bilateral)
%   'interneuron' - Task-modulated but fails classic criteria
%   'excluded'    - Below thresholds or not modulated
%
% For SNc sessions, call screen_snc_neurons for arrayROC-based classification:
%   'modulated'   - FR > 1 sp/s AND not sparse AND task modulated (3 consecutive sig. bins)
%   'excluded'    - Fails any of the three criteria

if strcmp(brain_area, 'SC') && ~isempty(project_root)
    fprintf('apply_neuron_screening: Running SC-specific classification...\n');

    try
        % Call screen_sc_neurons to get the detailed neuron_class.
        % Also capture updated session_data with cached_friedman_results
        % Request 6 outputs to get v2 classification. The first output
        % (sc_selected) is screen_sc_neurons' own inclusion decision.
        [sc_selected, ~, ~, sc_neuron_class, session_data, ...
            sc_neuron_class_v2] = ...
            screen_sc_neurons(session_data, ...
            project_root, condition_defs);

        % Transfer neuron_class (v1) to screening_info
        for i = 1:nNeurons
            if i <= length(sc_neuron_class)
                screening_info(i).neuron_class = ...
                    sc_neuron_class{i};
            end
        end

        % Adopt screen_sc_neurons' selection as the inclusion decision
        % for SC, so that a selected neuron always has a 'classic' or
        % 'interneuron' label and vice versa. Inclusion and
        % classification then share one definition of task modulation
        % (spatial Friedman/ranksum over visual/delay/saccade) rather
        % than the earlier split where inclusion used the baseline-
        % deviation test in compute_task_modulation.
        selected_neurons = sc_selected(:);
        for i = 1:nNeurons
            if i <= length(sc_selected)
                screening_info(i).included = sc_selected(i);
                if ~sc_selected(i)
                    screening_info(i).exclusion_reason = ...
                        'Not classic/interneuron (screen_sc)';
                else
                    screening_info(i).exclusion_reason = '';
                end
            end
        end

        % Transfer neuron_class_v2 to screening_info
        if exist('sc_neuron_class_v2', 'var') ...
                && ~isempty(sc_neuron_class_v2)
            for i = 1:nNeurons
                if i <= length(sc_neuron_class_v2)
                    screening_info(i).neuron_class_v2 ...
                        = sc_neuron_class_v2{i};
                end
            end
        end

        % Count and report SC classifications
        n_classic = sum(strcmp( ...
            sc_neuron_class, 'classic'));
        n_interneuron = sum(strcmp( ...
            sc_neuron_class, 'interneuron'));
        n_sc_excluded = sum(strcmp( ...
            sc_neuron_class, 'excluded'));

        fprintf(['apply_neuron_screening: SC ' ...
            'classification: %d classic, ' ...
            '%d interneuron, %d excluded\n'], ...
            n_classic, n_interneuron, ...
            n_sc_excluded);

    catch ME
        warning('apply_neuron_screening: SC classification failed: %s', ME.message);
        % Fall back to binary classification for SC
        for i = 1:nNeurons
            if screening_info(i).included
                screening_info(i).neuron_class = 'modulated';
            else
                screening_info(i).neuron_class = 'excluded';
            end
        end
    end

elseif strcmp(brain_area, 'SC') && isempty(project_root)
    % SC session but no project_root - use binary classification
    fprintf('apply_neuron_screening: SC session but no project_root provided. Using binary classification.\n');
    for i = 1:nNeurons
        if screening_info(i).included
            screening_info(i).neuron_class = 'modulated';
        else
            screening_info(i).neuron_class = 'excluded';
        end
    end

elseif strcmp(brain_area, 'SNc') && ~isempty(project_root)
    % SNc session - use screen_snc_neurons for arrayROC-based classification
    fprintf('apply_neuron_screening: Running SNc-specific classification...\n');

    try
        % Call screen_snc_neurons to get the detailed classification
        % This function implements the three-part test with arrayROC-based modulation
        % Pass modulation_info to avoid recomputing (significant time savings)
        [snc_selected, snc_screening_info, snc_neuron_class] = ...
            screen_snc_neurons(session_data, core_data, condition_defs, project_root, modulation_info);

        % Transfer neuron_class and update screening_info with SNc-specific results
        for i = 1:nNeurons
            if i <= length(snc_neuron_class)
                screening_info(i).neuron_class = snc_neuron_class{i};
                % Update included status based on screen_snc_neurons result
                screening_info(i).included = snc_selected(i);
                % Copy any additional modulation details from SNc screening
                if isfield(snc_screening_info(i), 'modulated_at')
                    screening_info(i).modulated_at = snc_screening_info(i).modulated_at;
                end
                if isfield(snc_screening_info(i), 'exclusion_reason') && ...
                        ~isempty(snc_screening_info(i).exclusion_reason)
                    screening_info(i).exclusion_reason = snc_screening_info(i).exclusion_reason;
                end
            end
        end

        % Update selected_neurons to match screen_snc_neurons result
        selected_neurons = snc_selected;

        % Count and report SNc classifications
        n_modulated = sum(strcmp(snc_neuron_class, 'modulated'));
        n_snc_excluded = sum(strcmp(snc_neuron_class, 'excluded'));

        fprintf('apply_neuron_screening: SNc classification: %d modulated, %d excluded\n', ...
            n_modulated, n_snc_excluded);

    catch ME
        warning('apply_neuron_screening: SNc classification failed: %s', ME.message);
        fprintf('Error details: %s\n', ME.message);
        % Fall back to binary classification for SNc
        for i = 1:nNeurons
            if screening_info(i).included
                screening_info(i).neuron_class = 'modulated';
            else
                screening_info(i).neuron_class = 'excluded';
            end
        end
    end

else
    % SNc without project_root or other brain areas - use binary classification
    for i = 1:nNeurons
        if screening_info(i).included
            screening_info(i).neuron_class = 'modulated';
        else
            screening_info(i).neuron_class = 'excluded';
        end
    end
end

%% Stability Gate
% Exclude neurons whose stable period covers too little of the in-scope
% trials (e.g. units isolated only during other tasks, or lost early in
% the gSac_4factors block). Applied after area-specific classification
% so that it holds for both SC and SNc.
for i = 1:nNeurons
    screening_info(i).stability_coverage = stability_info(i).coverage;
    screening_info(i).passes_stability = ...
        stability_info(i).is_stable && ...
        stability_info(i).coverage >= stability_params.min_coverage;
end
if use_strict_screening && stability_params.min_coverage > 0
    n_gated = 0;
    for i = 1:nNeurons
        if screening_info(i).included && ~screening_info(i).passes_stability
            screening_info(i).included = false;
            screening_info(i).neuron_class = 'excluded';
            screening_info(i).exclusion_reason = sprintf(...
                'Unstable in %s (coverage %.0f%%)', ...
                stability_params.trial_scope, ...
                100 * screening_info(i).stability_coverage);
            n_gated = n_gated + 1;
        end
    end
    selected_neurons = [screening_info.included]';
    fprintf(['apply_neuron_screening: Stability gate (coverage >= ' ...
        '%.0f%% of %s trials) excluded %d otherwise-included ' ...
        'neurons. Final: %d/%d included.\n'], ...
        100 * stability_params.min_coverage, ...
        stability_params.trial_scope, n_gated, ...
        sum(selected_neurons), nNeurons);
end

end

function scope_trials = get_stability_scope_trials(session_data, scope)
% GET_STABILITY_SCOPE_TRIALS Session trial indices used for stability
% detection: successful gSac_4factors trials, or all trials.
nTrials = numel(session_data.eventTimes.trialBegin);
switch scope
    case 'gSac_4factors'
        codes = initCodes();
        ti = session_data.trialInfo;
        is_4f = ti.taskCode(:) == codes.uniqueTaskCode_gSac_4factors;
        is_success = ti.trialEndState(:) == 21 & ~cellfun(@isempty, ...
            session_data.eventTimes.rewardCell(:));
        scope_trials = find(is_4f & is_success);
    case 'all'
        scope_trials = (1:nTrials)';
    otherwise
        error('apply_neuron_screening:badScope', ...
            'Unknown stability trial_scope "%s".', scope);
end
end
