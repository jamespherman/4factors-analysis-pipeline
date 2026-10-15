function aggregate_data = aggregate_screening_data(manifest, project_root)
% AGGREGATE_SCREENING_DATA Collect screening metrics from all sessions.
%
%   aggregate_data = AGGREGATE_SCREENING_DATA(manifest, project_root)
%
%   Loads screening_info from all processed sessions and combines into
%   a single structure for aggregate visualization.
%
%   INPUTS:
%       manifest     - Table from session_manifest.csv
%       project_root - Root directory of the project
%
%   OUTPUT:
%       aggregate_data - Struct containing all metrics across sessions:
%           .session_ids       - Cell array of session IDs (per neuron)
%           .brain_areas       - Cell array ('SC' or 'SNc')
%           .subjects          - Cell array of subject names
%           .cluster_ids       - Cluster IDs within session
%           .mean_firing_rates - Mean FR (sp/s)
%           .sparsity          - Proportion empty bins
%           .stability_scores  - Stability scores (0-1)
%           .presence_ratios   - Presence ratios (0-1)
%           .stable_period_frac- Fraction of trials in stable period
%           .is_modulated      - Logical: shows task modulation
%           .is_included       - Logical: passed screening
%           .exclusion_reasons - Cell array of strings
%           .neuron_classes    - Cell array of class labels
%           .fr_thresholds     - FR threshold used for each neuron
%           .session_indices   - Index into unique session list
%           .unique_sessions   - Cell array of unique session IDs
%           .session_summary   - Table with per-session summary stats
%
%   Author: Claude Code
%   Date: 2026-01-23

fprintf('aggregate_screening_data: Collecting screening data from %d sessions...\n', height(manifest));

% Initialize storage
all_session_ids = {};
all_brain_areas = {};
all_subjects = {};
all_cluster_ids = [];
all_mean_frs = [];
all_sparsity = [];
all_stability_scores = [];
all_presence_ratios = [];
all_stable_period_frac = [];
all_is_modulated = [];
all_is_included = [];
all_exclusion_reasons = {};
all_neuron_classes = {};
all_fr_thresholds = [];
all_session_indices = [];

% Per-session summary
session_summary = table();

unique_sessions = manifest.unique_id;
n_sessions = height(manifest);

for i = 1:n_sessions
    session_id = manifest.unique_id{i};
    brain_area = manifest.brain_area{i};
    subject = manifest.monkey{i};

    % Load session data
    local_processed_dir = fullfile(project_root, 'data', 'processed', session_id);
    local_session_data_path = fullfile(local_processed_dir, [session_id, '_session_data.mat']);

    if ~exist(local_session_data_path, 'file')
        fprintf('  Skipping %s: session data not found\n', session_id);
        continue;
    end

    fprintf('  Loading %s...\n', session_id);
    loaded = load(local_session_data_path, 'session_data');
    session_data = loaded.session_data;

    % Check for screening_info
    if ~isfield(session_data, 'analysis') || ~isfield(session_data.analysis, 'screening_info')
        fprintf('  Skipping %s: no screening_info found\n', session_id);
        continue;
    end

    screening_info = session_data.analysis.screening_info;
    nNeurons = length(screening_info);

    % Determine FR threshold
    if strcmpi(brain_area, 'SNc')
        fr_threshold = 1.0;
    else
        fr_threshold = 3.5;
    end

    % Extract stability info if available
    has_stability = isfield(session_data, 'stability_info') && ~isempty(session_data.stability_info);
    if has_stability
        stability_info = session_data.stability_info;
        nTrials = length(session_data.eventTimes.trialBegin);
    end

    % Extract metrics for each neuron
    for j = 1:nNeurons
        all_session_ids{end+1} = session_id;
        all_brain_areas{end+1} = brain_area;
        all_subjects{end+1} = subject;
        all_session_indices(end+1) = i;

        % Basic metrics
        all_cluster_ids(end+1) = screening_info(j).cluster_id;
        all_mean_frs(end+1) = screening_info(j).mean_firing_rate;
        all_sparsity(end+1) = screening_info(j).proportion_empty_bins;
        all_is_included(end+1) = screening_info(j).included;
        all_fr_thresholds(end+1) = fr_threshold;

        % Modulation
        if isfield(screening_info(j), 'is_modulated')
            all_is_modulated(end+1) = screening_info(j).is_modulated;
        else
            all_is_modulated(end+1) = false;
        end

        % Exclusion reason
        if isfield(screening_info(j), 'exclusion_reason') && ~isempty(screening_info(j).exclusion_reason)
            all_exclusion_reasons{end+1} = screening_info(j).exclusion_reason;
        else
            all_exclusion_reasons{end+1} = '';
        end

        % Neuron class
        if isfield(screening_info(j), 'neuron_class')
            all_neuron_classes{end+1} = screening_info(j).neuron_class;
        else
            all_neuron_classes{end+1} = '';
        end

        % Stability metrics
        if has_stability && j <= length(stability_info)
            all_stability_scores(end+1) = stability_info(j).stability_score;
            all_presence_ratios(end+1) = stability_info(j).presence_ratio;
            all_stable_period_frac(end+1) = stability_info(j).n_stable_trials / nTrials;
        else
            all_stability_scores(end+1) = NaN;
            all_presence_ratios(end+1) = NaN;
            all_stable_period_frac(end+1) = NaN;
        end
    end

    % Compute session summary
    session_row = table();
    session_row.session_id = {session_id};
    session_row.brain_area = {brain_area};
    session_row.subject = {subject};
    session_row.n_neurons = nNeurons;
    session_row.n_included = sum([screening_info.included]);
    session_row.pct_included = 100 * session_row.n_included / nNeurons;
    session_row.mean_fr = mean([screening_info.mean_firing_rate]);
    session_row.median_fr = median([screening_info.mean_firing_rate]);
    session_row.mean_sparsity = mean([screening_info.proportion_empty_bins]);
    session_row.pct_pass_fr = 100 * sum([screening_info.passes_fr_threshold]) / nNeurons;
    session_row.pct_pass_sparsity = 100 * sum([screening_info.passes_sparsity_threshold]) / nNeurons;

    % Modulation rate
    n_modulated = 0;
    for j = 1:nNeurons
        if isfield(screening_info(j), 'is_modulated') && screening_info(j).is_modulated
            n_modulated = n_modulated + 1;
        end
    end
    session_row.pct_modulated = 100 * n_modulated / nNeurons;

    % Stability summary
    if has_stability
        session_row.mean_stability_score = mean([stability_info.stability_score]);
        session_row.pct_stable = 100 * sum([stability_info.is_stable]) / nNeurons;
    else
        session_row.mean_stability_score = NaN;
        session_row.pct_stable = NaN;
    end

    session_summary = [session_summary; session_row];
end

% Convert to column vectors
all_mean_frs = all_mean_frs(:);
all_sparsity = all_sparsity(:);
all_stability_scores = all_stability_scores(:);
all_presence_ratios = all_presence_ratios(:);
all_stable_period_frac = all_stable_period_frac(:);
all_is_modulated = all_is_modulated(:);
all_is_included = all_is_included(:);
all_cluster_ids = all_cluster_ids(:);
all_fr_thresholds = all_fr_thresholds(:);
all_session_indices = all_session_indices(:);

% Build output structure
aggregate_data = struct();
aggregate_data.session_ids = all_session_ids(:);
aggregate_data.brain_areas = all_brain_areas(:);
aggregate_data.subjects = all_subjects(:);
aggregate_data.cluster_ids = all_cluster_ids;
aggregate_data.mean_firing_rates = all_mean_frs;
aggregate_data.sparsity = all_sparsity;
aggregate_data.stability_scores = all_stability_scores;
aggregate_data.presence_ratios = all_presence_ratios;
aggregate_data.stable_period_frac = all_stable_period_frac;
aggregate_data.is_modulated = logical(all_is_modulated);
aggregate_data.is_included = logical(all_is_included);
aggregate_data.exclusion_reasons = all_exclusion_reasons(:);
aggregate_data.neuron_classes = all_neuron_classes(:);
aggregate_data.fr_thresholds = all_fr_thresholds;
aggregate_data.session_indices = all_session_indices;
aggregate_data.unique_sessions = unique_sessions;
aggregate_data.session_summary = session_summary;

fprintf('aggregate_screening_data: Collected %d neurons from %d sessions.\n', ...
    length(all_mean_frs), height(session_summary));
fprintf('  Total included: %d (%.1f%%)\n', sum(all_is_included), 100*mean(all_is_included));

end
