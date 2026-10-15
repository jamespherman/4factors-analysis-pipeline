function clustering = cluster_sc_neurons( ...
    feature_matrix, feature_names, cluster_ids, opts)
% cluster_sc_neurons - Unsupervised clustering of SC
% neurons based on response features.
%
% Determines optimal number of clusters using silhouette
% analysis and gap statistic, then assigns neurons via
% k-means. Also runs hierarchical clustering for
% comparison and bootstrap for stability assessment.
%
% INPUTS:
%   feature_matrix - [N x F] z-scored feature matrix
%                    (only quality-passing neurons)
%   feature_names  - {1 x F} cell of feature names
%   cluster_ids    - [N x 1] neuron cluster IDs
%   opts           - (Optional) struct with fields:
%       .max_k         Max clusters to test (default 10)
%       .min_per_cluster  Min neurons per cluster (def 15)
%       .n_bootstrap   Number of bootstrap iterations
%                      (default 100)
%       .random_seed   RNG seed (default 42)
%       .force_k       Force a specific k (default [])
%
% OUTPUT:
%   clustering - Struct with fields:
%       .assignments     [N x 1] cluster index (1..k)
%       .k               optimal k
%       .centroids       [k x F] cluster centroids
%       .silhouette_avg  mean silhouette score
%       .silhouette_vals [N x 1] per-neuron silhouette
%       .eval_silhouette evalclusters result (silhouette)
%       .eval_gap        evalclusters result (gap)
%       .bic_values      [1 x max_k] BIC for GMM
%       .pca_scores      [N x nPC] PCA-transformed
%       .pca_explained   variance explained per PC
%       .pca_coeff       [F x nPC] PCA loadings
%       .linkage_tree    hierarchical linkage matrix
%       .hier_assignments [N x 1] hierarchical result
%       .bootstrap_stability [N x 1] proportion of
%                        bootstrap runs agreeing
%       .labels          {k x 1} cluster labels
%       .feature_names   {1 x F} passed through
%       .cluster_ids     [N x 1] passed through
%
% Author: Claude Code
% Date: 2026-02-06

%% Parse options
if nargin < 4 || isempty(opts)
    opts = struct();
end
max_k = getfield_default(opts, 'max_k', 10);
min_per_cluster = getfield_default( ...
    opts, 'min_per_cluster', 10);
n_bootstrap = getfield_default( ...
    opts, 'n_bootstrap', 100);
random_seed = getfield_default( ...
    opts, 'random_seed', 42);
force_k = getfield_default(opts, 'force_k', []);

[N, F] = size(feature_matrix);
fprintf(['cluster_sc_neurons: %d neurons, ' ...
    '%d features.\n'], N, F);

if N < 30
    warning(['cluster_sc_neurons: Only %d neurons. ' ...
        'Results may be unreliable.'], N);
end

%% Set random seed for reproducibility
rng(random_seed);

%% Step 1: PCA for visualization
% (also used for clustering if many correlated features)
[pca_coeff, pca_scores, ~, ~, pca_explained] = ...
    pca(feature_matrix);

% Determine number of PCs capturing 90% variance
cumvar = cumsum(pca_explained);
n_pc = find(cumvar >= 90, 1);
if isempty(n_pc)
    n_pc = F;
end
fprintf(['cluster_sc_neurons: %d PCs capture ' ...
    '%.1f%% variance.\n'], n_pc, cumvar(n_pc));

% Check for high correlations among features
corr_matrix = corrcoef(feature_matrix);
% Count pairs with |r| > 0.8 (excluding diagonal)
high_corr = sum(abs(corr_matrix(:)) > 0.8) ...
    - F;  % subtract diagonal
high_corr_pairs = high_corr / 2;

% Use PCA-reduced features for clustering if many
% correlated features
if high_corr_pairs > F
    fprintf(['cluster_sc_neurons: %d highly ' ...
        'correlated pairs. Using PCA features ' ...
        'for clustering.\n'], high_corr_pairs);
    clust_features = pca_scores(:, 1:n_pc);
else
    clust_features = feature_matrix;
end

%% Step 2: Determine optimal k
% Limit max_k based on sample size
max_k = min(max_k, floor(N / min_per_cluster));
max_k = max(max_k, 2);  % at least test k=2

k_range = 2:max_k;
fprintf(['cluster_sc_neurons: Testing k = %d ' ...
    'to %d.\n'], k_range(1), k_range(end));

% Silhouette criterion
eval_sil = evalclusters( ...
    clust_features, 'kmeans', 'silhouette', ...
    'KList', k_range);

% Gap statistic
eval_gap = evalclusters( ...
    clust_features, 'kmeans', 'gap', ...
    'KList', k_range);

% BIC from Gaussian Mixture Models
bic_values = nan(1, max_k);
for k = k_range
    try
        gm = fitgmdist(clust_features, k, ...
            'RegularizationValue', 0.01, ...
            'Replicates', 5, ...
            'Options', statset('MaxIter', 200));
        bic_values(k) = gm.BIC;
    catch
        % GMM can fail for small k or
        % ill-conditioned data
        bic_values(k) = Inf;
    end
end

% Find BIC-optimal k
valid_bic = find(isfinite(bic_values));
if ~isempty(valid_bic)
    [~, bic_opt_idx] = min( ...
        bic_values(valid_bic));
    bic_optimal_k = valid_bic(bic_opt_idx);
else
    bic_optimal_k = NaN;
end

% Choose optimal k
if ~isempty(force_k)
    optimal_k = force_k;
    fprintf(['cluster_sc_neurons: Using forced ' ...
        'k = %d.\n'], optimal_k);
else
    % Prefer silhouette, but cross-check
    optimal_k = eval_sil.OptimalK;
    fprintf(['cluster_sc_neurons: Silhouette ' ...
        'optimal k=%d, Gap optimal k=%d, ' ...
        'BIC optimal k=%d.\n'], ...
        eval_sil.OptimalK, ...
        eval_gap.OptimalK, ...
        bic_optimal_k);

    % If all three agree, great. Otherwise use
    % silhouette as primary criterion.
    if eval_sil.OptimalK == eval_gap.OptimalK
        optimal_k = eval_sil.OptimalK;
    end
end

fprintf(['cluster_sc_neurons: Selected k = %d.\n'], ...
    optimal_k);

%% Step 3: Final k-means clustering
[assignments, centroids, sumd] = kmeans( ...
    clust_features, optimal_k, ...
    'Replicates', 20, ...
    'MaxIter', 500, ...
    'Distance', 'sqeuclidean', ...
    'Display', 'off');

%% Step 3b: Enforce minimum cluster size
% Merge small clusters into nearest neighbor
merge_happened = true;
while merge_happened
    merge_happened = false;
    cluster_counts = histcounts( ...
        assignments, 1:(optimal_k+1));
    small_clusters = find( ...
        cluster_counts < min_per_cluster);
    if isempty(small_clusters)
        break;
    end
    % Merge the smallest cluster
    [~, sm_idx] = min( ...
        cluster_counts(small_clusters));
    sm_c = small_clusters(sm_idx);
    fprintf(['cluster_sc_neurons: Cluster %d ' ...
        'has %d neurons (< %d). Merging.\n'], ...
        sm_c, cluster_counts(sm_c), ...
        min_per_cluster);

    % Find nearest centroid to merge into
    dists = pdist2( ...
        centroids(sm_c, :), centroids);
    dists(sm_c) = Inf;
    [~, merge_target] = min(dists);
    fprintf(['  -> Merging into cluster %d.\n'], ...
        merge_target);

    % Reassign neurons
    assignments(assignments == sm_c) = ...
        merge_target;

    % Renumber clusters to be contiguous
    old_k = optimal_k;
    optimal_k = optimal_k - 1;
    new_assignments = zeros(size(assignments));
    remaining = setdiff(1:old_k, sm_c);
    for i_r = 1:length(remaining)
        new_assignments( ...
            assignments == remaining(i_r)) = i_r;
    end
    assignments = new_assignments;

    % Recompute centroids
    centroids = zeros(optimal_k, ...
        size(clust_features, 2));
    for c = 1:optimal_k
        centroids(c, :) = mean( ...
            clust_features( ...
            assignments == c, :), 1);
    end

    merge_happened = true;
end

fprintf(['cluster_sc_neurons: Final k = %d ' ...
    'after merging.\n'], optimal_k);

%% Step 4: Hierarchical clustering for comparison
linkage_tree = linkage(clust_features, 'ward');
hier_assignments = cluster(linkage_tree, ...
    'maxclust', optimal_k);

%% Step 5: Silhouette scores
[sil_vals, ~] = silhouette( ...
    clust_features, assignments);
sil_avg = mean(sil_vals);
fprintf(['cluster_sc_neurons: Mean silhouette = ' ...
    '%.3f.\n'], sil_avg);

% Per-cluster silhouette
for k = 1:optimal_k
    k_mask = assignments == k;
    k_sil = mean(sil_vals(k_mask));
    fprintf(['  Cluster %d: %d neurons, ' ...
        'silhouette = %.3f\n'], ...
        k, sum(k_mask), k_sil);
end

%% Step 6: Bootstrap stability
fprintf(['cluster_sc_neurons: Running %d ' ...
    'bootstrap iterations...\n'], n_bootstrap);
boot_agree = zeros(N, 1);
boot_appearances = zeros(N, 1);

for b = 1:n_bootstrap
    % Resample with replacement
    boot_idx = randsample(N, N, true);
    boot_features = clust_features(boot_idx, :);

    % Track which neurons appeared
    in_bag = unique(boot_idx);
    boot_appearances(in_bag) = ...
        boot_appearances(in_bag) + 1;

    % Cluster the bootstrap sample
    try
        boot_assign = kmeans( ...
            boot_features, optimal_k, ...
            'Start', centroids, ...
            'MaxIter', 100, ...
            'Display', 'off');
    catch
        continue;
    end

    % Map bootstrap clusters to original using
    % contingency table (Hungarian-like greedy)
    mapping = match_clusters( ...
        assignments(boot_idx), boot_assign, ...
        optimal_k);

    % Check agreement for in-bag neurons
    mapped_boot = mapping(boot_assign);
    for i_n = 1:length(in_bag)
        neuron = in_bag(i_n);
        matches = find(boot_idx == neuron);
        if all(mapped_boot(matches) == ...
                assignments(neuron))
            boot_agree(neuron) = ...
                boot_agree(neuron) + 1;
        end
    end
end

% Normalize by per-neuron appearance count
bootstrap_stability = boot_agree ./ ...
    max(boot_appearances, 1);

fprintf(['cluster_sc_neurons: Mean bootstrap ' ...
    'stability = %.1f%%.\n'], ...
    mean(bootstrap_stability) * 100);

%% Step 7: Initialize labels
labels = cell(optimal_k, 1);
for k = 1:optimal_k
    labels{k} = sprintf('cluster_%d', k);
end

%% Pack output
clustering = struct();
clustering.assignments = assignments;
clustering.k = optimal_k;
clustering.centroids = centroids;
clustering.silhouette_avg = sil_avg;
clustering.silhouette_vals = sil_vals;
clustering.eval_silhouette = eval_sil;
clustering.eval_gap = eval_gap;
clustering.bic_values = bic_values;
clustering.pca_scores = pca_scores;
clustering.pca_explained = pca_explained;
clustering.pca_coeff = pca_coeff;
clustering.linkage_tree = linkage_tree;
clustering.hier_assignments = hier_assignments;
clustering.bootstrap_stability = bootstrap_stability;
clustering.labels = labels;
clustering.feature_names = feature_names;
clustering.cluster_ids = cluster_ids;
clustering.clust_features = clust_features;
clustering.sumd = sumd;

end

%% Helper: match clusters between two assignments
function mapping = match_clusters( ...
    ref_assign, test_assign, k)
% Greedy matching: for each test cluster, find the
% reference cluster with highest overlap.
mapping = zeros(k, 1);
used_ref = false(k, 1);

% Build contingency table
contingency = zeros(k, k);
for i = 1:k
    for j = 1:k
        contingency(i, j) = sum( ...
            ref_assign == i & test_assign == j);
    end
end

% Greedy assignment (largest overlap first)
for iter = 1:k
    [~, max_idx] = max(contingency(:));
    [ref_k, test_k] = ind2sub([k, k], max_idx);
    mapping(test_k) = ref_k;
    used_ref(ref_k) = true;
    % Zero out this row and column
    contingency(ref_k, :) = 0;
    contingency(:, test_k) = 0;
end

% Handle unmapped clusters (if any)
unmapped = find(mapping == 0);
unused = find(~used_ref);
for i = 1:length(unmapped)
    if i <= length(unused)
        mapping(unmapped(i)) = unused(i);
    else
        mapping(unmapped(i)) = 1;
    end
end
end

%% Helper: get struct field with default
function val = getfield_default(s, field, default)
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end
