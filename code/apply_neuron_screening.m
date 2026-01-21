function [selected_neurons, screening_info] = apply_neuron_screening(session_data, core_data, condition_defs, project_root)
% APPLY_NEURON_SCREENING Apply refined neuron screening criteria.
%
%   [selected_neurons, screening_info] = APPLY_NEURON_SCREENING(session_data,
%       core_data, condition_defs)
%   [selected_neurons, screening_info] = APPLY_NEURON_SCREENING(session_data,
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
end

%% Criterion 1 & 3: Mean Firing Rate and Sparsity (computed over ACTIVE PERIOD)
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
        continue;
    end

    % Define active period for this neuron (first spike to last spike)
    active_start = spike_times(1);
    active_end = spike_times(end);
    active_duration = active_end - active_start;

    % Store active period info
    screening_info(i).active_start = active_start;
    screening_info(i).active_end = active_end;
    screening_info(i).active_duration = active_duration;

    % Handle edge case: very short active period
    if active_duration < 1  % Less than 1 second
        screening_info(i).mean_firing_rate = n_spikes;  % Treat as if 1 second
        screening_info(i).proportion_empty_bins = 0;    % Can't compute meaningfully
        screening_info(i).passes_fr_threshold = n_spikes >= fr_threshold;
        screening_info(i).passes_sparsity_threshold = true;  % Benefit of doubt
        continue;
    end

    % Compute mean firing rate over active period only
    mean_fr = n_spikes / active_duration;
    screening_info(i).mean_firing_rate = mean_fr;
    screening_info(i).passes_fr_threshold = mean_fr >= fr_threshold;

    % Compute sparsity over active period only (100ms bins)
    edges = active_start : 0.1 : active_end;
    if length(edges) < 2
        % Not enough bins
        screening_info(i).proportion_empty_bins = 0;
        screening_info(i).passes_sparsity_threshold = true;
    else
        counts = histcounts(spike_times, edges);
        proportion_empty = sum(counts == 0) / numel(counts);
        screening_info(i).proportion_empty_bins = proportion_empty;
        screening_info(i).passes_sparsity_threshold = proportion_empty < max_proportion_empty;
    end
end

n_pass_fr = sum([screening_info.passes_fr_threshold]);
n_pass_sparsity = sum([screening_info.passes_sparsity_threshold]);
fprintf('apply_neuron_screening: %d/%d neurons (%.1f%%) pass FR threshold (mean FR >= %.1f sp/s over active period).\n', ...
    n_pass_fr, nNeurons, 100 * n_pass_fr / nNeurons, fr_threshold);
fprintf('apply_neuron_screening: %d/%d neurons (%.1f%%) pass sparsity threshold (<%.0f%% empty bins over active period).\n', ...
    n_pass_sparsity, nNeurons, 100 * n_pass_sparsity / nNeurons, max_proportion_empty * 100);

%% Criterion 2: Task Modulation
% Call compute_task_modulation to check for modulation
modulation_info = compute_task_modulation(session_data, core_data, condition_defs);

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
        % Call screen_sc_neurons to get the detailed neuron_class
        [~, ~, ~, sc_neuron_class] = screen_sc_neurons(session_data, project_root, condition_defs);

        % Transfer neuron_class to screening_info
        for i = 1:nNeurons
            if i <= length(sc_neuron_class)
                screening_info(i).neuron_class = sc_neuron_class{i};
            end
        end

        % Count and report SC classifications
        n_classic = sum(strcmp(sc_neuron_class, 'classic'));
        n_interneuron = sum(strcmp(sc_neuron_class, 'interneuron'));
        n_sc_excluded = sum(strcmp(sc_neuron_class, 'excluded'));

        fprintf('apply_neuron_screening: SC classification: %d classic, %d interneuron, %d excluded\n', ...
            n_classic, n_interneuron, n_sc_excluded);

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
        [snc_selected, snc_screening_info, snc_neuron_class] = ...
            screen_snc_neurons(session_data, core_data, condition_defs, project_root);

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

end
