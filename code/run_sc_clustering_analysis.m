%% run_sc_clustering_analysis.m
%
% Standalone driver for SC neuron clustering (v2
% classification). Loads all SC sessions, extracts
% features, runs clustering, generates visualizations,
% and validates.
%
% This script is NOT part of the main pipeline. It is
% used for development and validation of the data-driven
% SC neuron classification.
%
% Usage:
%   matlab -batch "run('code/run_sc_clustering_analysis.m')"
%
% Author: Claude Code
% Date: 2026-02-06

clear; clc; close all;

%% Setup
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

tic;
giveFeed = @(x) disp( ...
    [num2str(round(toc, 1)) 's - ' x]);

%% Load manifest and identify SC sessions
manifest_path = fullfile(project_root, ...
    'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

sc_mask = strcmp(manifest.brain_area, 'SC');
sc_manifest = manifest(sc_mask, :);
n_sessions = height(sc_manifest);

giveFeed(sprintf('Found %d SC sessions.', ...
    n_sessions));

%% Phase 1: Extract features from each session
all_features_raw = [];
all_features_z = [];
all_cluster_ids = [];
all_session_ids = {};
all_v1_classes = {};
all_passes_quality = logical([]);

for i_s = 1:n_sessions
    session_id = sc_manifest.unique_id{i_s};
    giveFeed(sprintf( ...
        'Processing session %d/%d: %s', ...
        i_s, n_sessions, session_id));

    %% Load session data
    local_path = fullfile(project_root, ...
        'data', 'processed', session_id, ...
        [session_id '_session_data.mat']);

    if ~exist(local_path, 'file')
        % Try OneDrive
        one_drive_path = findOneDrive;
        local_path = fullfile(one_drive_path, ...
            'Neuronal Data Analysis', ...
            session_id, ...
            [session_id '_session_data.mat']);
    end

    if ~exist(local_path, 'file')
        warning('Session data not found: %s', ...
            session_id);
        continue;
    end

    load(local_path, 'session_data');

    %% Check for cached_friedman_results
    if ~isfield(session_data, ...
            'cached_friedman_results') || ...
            isempty( ...
            session_data.cached_friedman_results)
        warning(['No cached_friedman_results ' ...
            'for %s. Run main pipeline first.'], ...
            session_id);
        continue;
    end

    %% Get scSide (need to re-derive or use
    %  analysis results)
    scSide = 'unknown';
    if isfield(session_data, 'analysis') && ...
            isfield(session_data.analysis, ...
            'screening_info')
        % Look for scSide in the screening output
        % (not directly stored, derive from targets)
    end

    % Derive scSide from gSac_jph target locations
    codes = initCodes();
    trialInfo = session_data.trialInfo;
    eventTimes = session_data.eventTimes;

    jph_trials = find( ...
        trialInfo.taskCode == ...
        codes.uniqueTaskCode_gSac_jph & ...
        eventTimes.targetReillum > 0 & ...
        eventTimes.pdsReward > 0);

    if length(jph_trials) > 5
        thetas = trialInfo.targetTheta( ...
            jph_trials) / 10;
        left_vf = sum( ...
            thetas > 90 & thetas < 270);
        right_vf = sum( ...
            thetas < 90 | thetas > 270);
        if left_vf > right_vf
            scSide = 'right';
        else
            scSide = 'left';
        end
    else
        % Fallback: grid hole
        if isfield(session_data, 'metadata') ...
                && isfield( ...
                session_data.metadata, ...
                'grid_hole')
            gh = session_data.metadata.grid_hole;
            coords = sscanf(gh, '(%f, %f)');
            if ~isempty(coords)
                if coords(1) < 0
                    scSide = 'left';
                elseif coords(1) > 0
                    scSide = 'right';
                end
            end
        end
    end

    %% Extract features
    features = extract_sc_response_features( ...
        session_data, scSide);

    %% Get v1 classification
    v1_class = repmat({'excluded'}, ...
        length(features.cluster_ids), 1);
    if isfield(session_data, 'analysis') && ...
            isfield(session_data.analysis, ...
            'screening_info')
        si = session_data.analysis.screening_info;
        for i_n = 1:min(length(si), ...
                length(features.cluster_ids))
            if isfield(si(i_n), 'neuron_class') ...
                    && ~isempty( ...
                    si(i_n).neuron_class)
                v1_class{i_n} = ...
                    si(i_n).neuron_class;
            end
        end
    end

    %% Accumulate across sessions
    n_neurons = length(features.cluster_ids);
    all_features_raw = [all_features_raw; ...
        features.raw_features]; %#ok<AGROW>
    all_cluster_ids = [all_cluster_ids; ...
        features.cluster_ids]; %#ok<AGROW>
    all_passes_quality = [all_passes_quality; ...
        features.passes_quality]; %#ok<AGROW>
    all_v1_classes = [all_v1_classes; ...
        v1_class]; %#ok<AGROW>

    session_labels = repmat({session_id}, ...
        n_neurons, 1);
    all_session_ids = [all_session_ids; ...
        session_labels]; %#ok<AGROW>
end

giveFeed(sprintf( ...
    'Feature extraction complete. %d total neurons.', ...
    length(all_cluster_ids)));

%% Filter to quality-passing neurons only
quality_mask = all_passes_quality;
n_quality = sum(quality_mask);
giveFeed(sprintf( ...
    '%d / %d neurons pass quality gates.', ...
    n_quality, length(all_cluster_ids)));

if n_quality < 30
    error(['Too few quality-passing neurons ' ...
        '(%d). Need at least 30.'], n_quality);
end

% Extract quality-passing subset
q_features_raw = all_features_raw( ...
    quality_mask, :);
q_cluster_ids = all_cluster_ids(quality_mask);
q_session_ids = all_session_ids(quality_mask);
q_v1_classes = all_v1_classes(quality_mask);

% Feature names matching
% extract_sc_response_features output
feat_names = { ...
    'best_vis', '2nd_vis', ...
    '3rd_vis', 'worst_vis', ...
    'best_sac', '2nd_sac', ...
    '3rd_sac', 'worst_sac', ...
    'best_delay', 'worst_delay', ...
    'n_loc_sig_vis', 'n_loc_sig_sac', ...
    'n_loc_sig_delay', 'n_loc_suppressed', ...
    'vis_to_sac_ratio', 'delay_index', ...
    'has_suppression', 'has_mixed_polarity', ...
    'tuning_breadth'};
nFeatures = length(feat_names);

% Z-score normalize across quality neurons

q_features_z = q_features_raw;
for i_f = 1:nFeatures
    vals = q_features_raw(:, i_f);
    mu = mean(vals, 'omitnan');
    sd = std(vals, 'omitnan');
    if sd > 0
        q_features_z(:, i_f) = (vals - mu) / sd;
    else
        q_features_z(:, i_f) = 0;
    end
end

%% Phase 2: Clustering
giveFeed('Running clustering analysis...');
clustering = cluster_sc_neurons( ...
    q_features_z, feat_names, q_cluster_ids);

%% Phase 3: Visualization
giveFeed('Generating visualizations...');
features_agg = struct( ...
    'feature_matrix', q_features_z, ...
    'raw_features', q_features_raw, ...
    'feature_names', {feat_names});

visualize_sc_clusters( ...
    features_agg, clustering, ...
    q_v1_classes, project_root);

%% Phase 4: Validation
giveFeed('Running validation...');
validation = validate_sc_classification( ...
    features_agg, clustering, ...
    q_session_ids, q_v1_classes, project_root);

%% Print final summary
fprintf('\n');
fprintf('========================================\n');
fprintf('  SC Classification v2 Summary\n');
fprintf('========================================\n');
fprintf('Sessions:           %d\n', n_sessions);
fprintf('Total neurons:      %d\n', ...
    length(all_cluster_ids));
fprintf('Quality-passing:    %d\n', n_quality);
fprintf('Optimal k:          %d\n', clustering.k);
fprintf('Mean silhouette:    %.3f\n', ...
    clustering.silhouette_avg);
fprintf('Bootstrap stability: %.1f%%\n', ...
    mean(clustering.bootstrap_stability) * 100);
fprintf('LOO agreement:      %.1f%%\n', ...
    validation.loo_mean * 100);
fprintf('========================================\n');
fprintf('\nCluster assignments:\n');
for c = 1:clustering.k
    n_c = sum(clustering.assignments == c);
    fprintf('  %s: %d neurons (%.1f%%)\n', ...
        clustering.labels{c}, n_c, ...
        100 * n_c / n_quality);
end

%% Export CSV for expert validation
% Compute distance to centroid for each neuron
% (identifies exemplar neurons closest to center)
dist_to_centroid = nan(n_quality, 1);
for i_n = 1:n_quality
    c = clustering.assignments(i_n);
    dist_to_centroid(i_n) = sqrt(sum( ...
        (clustering.clust_features(i_n, :) - ...
        clustering.centroids(c, :)).^2));
end

% Build a table with all neuron info
neuron_table = table( ...
    q_session_ids, q_cluster_ids, ...
    clustering.assignments, ...
    dist_to_centroid, ...
    q_v1_classes, ...
    clustering.silhouette_vals, ...
    clustering.bootstrap_stability, ...
    'VariableNames', { ...
        'session_id', 'cluster_id', ...
        'v2_cluster', 'dist_to_centroid', ...
        'v1_class', 'silhouette', ...
        'boot_stability'});

% Add cluster label column
v2_label = cell(n_quality, 1);
for i_n = 1:n_quality
    v2_label{i_n} = clustering.labels{ ...
        clustering.assignments(i_n)};
end
neuron_table.v2_label = v2_label;

% Add PDF path column
pdf_paths = cell(n_quality, 1);
for i_n = 1:n_quality
    pdf_paths{i_n} = fullfile('figures', ...
        q_session_ids{i_n}, ...
        sprintf('%s_%03d_neuron_diagnostics.pdf', ...
        q_session_ids{i_n}, ...
        q_cluster_ids(i_n)));
end
neuron_table.pdf_path = pdf_paths;

% Sort by cluster, then distance to centroid
neuron_table = sortrows( ...
    neuron_table, {'v2_cluster', ...
    'dist_to_centroid'}, ...
    {'ascend', 'ascend'});

% Save CSV
csv_path = fullfile(project_root, 'figures', ...
    'sc_clustering_neuron_assignments.csv');
writetable(neuron_table, csv_path);
giveFeed(sprintf('Neuron assignments saved to %s', ...
    csv_path));

%% Print exemplar neurons for each cluster
% (5 closest to centroid per cluster)
n_exemplars = 5;
fprintf('\n=== Exemplar Neurons Per Cluster ===\n');
fprintf(['(closest to centroid - review their ' ...
    'diagnostic PDFs)\n\n']);

for c = 1:clustering.k
    c_rows = neuron_table( ...
        neuron_table.v2_cluster == c, :);
    n_show = min(n_exemplars, height(c_rows));

    fprintf('--- %s (n=%d) ---\n', ...
        clustering.labels{c}, height(c_rows));
    for i_n = 1:n_show
        fprintf('  %s cluster %d  (v1: %s, ', ...
            c_rows.session_id{i_n}, ...
            c_rows.cluster_id(i_n), ...
            c_rows.v1_class{i_n});
        fprintf('sil: %.2f)  %s\n', ...
            c_rows.silhouette(i_n), ...
            c_rows.pdf_path{i_n});
    end
    fprintf('\n');
end

%% Save results
output_path = fullfile(project_root, 'data', ...
    'sc_clustering_results.mat');
save(output_path, 'clustering', 'validation', ...
    'features_agg', 'q_cluster_ids', ...
    'q_session_ids', 'q_v1_classes', ...
    'feat_names', 'neuron_table', '-v7.3');
giveFeed(sprintf('Results saved to %s', ...
    output_path));

fprintf(['\n*** Next steps: ***\n']);
fprintf(['  1. Open sc_clustering_neuron_' ...
    'assignments.csv in figures/\n']);
fprintf(['  2. Review diagnostic PDFs for ' ...
    'exemplar neurons listed above\n']);
fprintf(['  3. Assign meaningful labels to ' ...
    'clusters in clustering.labels\n']);
fprintf(['  4. Re-save the .mat file with ' ...
    'updated labels\n']);
