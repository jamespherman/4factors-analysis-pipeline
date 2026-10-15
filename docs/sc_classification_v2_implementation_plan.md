# SC Neuron Classification v2 - Implementation Plan

**Created:** 2026-02-06
**Status:** Ready for implementation
**Prerequisite:** Open questions resolved (see Section 0)

---

## 0. Resolved Open Questions

These were discussed and resolved before implementation:

| # | Question | Decision |
|---|----------|----------|
| Q1 | Interpretability vs accuracy | Start with data-driven k, merge clusters if centroids not distinguishable and biology overlaps. "Classic" is ONE top-level category (V/VM/M/D are secondary attributes within it). Target 3-7 top-level categories. |
| Q2 | Treatment of excluded neurons | Do NOT cluster excluded neurons. Quality gates (FR>3, sparsity<=80%, task-modulated) are prerequisites. Only label exclusion reasons descriptively. |
| Q3 | Integration with existing analyses | Preserve v1 classification. Add v2 alongside. Both run in parallel. |
| Q4 | Manual overrides | Keep mechanism but don't invest time. If v2 fails smell test, user will build manual GUI. |
| Q5 | Version control | Phased: `neuron_class` (v1) preserved, `neuron_class_v2` added as new field. After validation, v2 becomes default. |

---

## 1. Architecture Overview

### Data Flow

```
cached_friedman_results (already in session_data)
        |
        v
extract_sc_response_features.m  [NEW - Phase 1]
        |
        v
feature_matrix (neurons x features)
        |
        v
cluster_sc_neurons.m  [NEW - Phase 2]
        |
        v
cluster_assignments + cluster_labels
        |
        v
screen_sc_neurons.m  [MODIFIED - Phase 5]
   (adds neuron_class_v2 field)
        |
        v
apply_neuron_screening.m  [MODIFIED - Phase 5]
   (propagates neuron_class_v2 to screening_info)
        |
        v
aggregate_analysis_results.m  [MODIFIED - Phase 5]
   (adds neuron_type_v2 column to neuron_metrics_table)
```

### Key Constraint
The clustering must be deterministic and reproducible. We will:
1. Use a fixed random seed for k-means
2. Store cluster centroids so new neurons can be assigned
3. Cache results in session_data (like cached_friedman_results)

---

## 2. Phase 1: Feature Extraction

### File: `code/extract_sc_response_features.m`

**Purpose:** Extract a feature vector for each SC neuron from
cached_friedman_results.

**Signature:**
```matlab
function features = extract_sc_response_features(session_data)
% Returns: features struct with fields:
%   .cluster_ids    [N x 1] cluster IDs
%   .feature_matrix [N x F] normalized feature values
%   .feature_names  {1 x F} cell of feature names
%   .raw_features   [N x F] before normalization
%   .passes_quality [N x 1] logical (FR, sparsity, modulated)
```

**Features to extract (per neuron):**

From `cached_friedman_results.group_mean_frs` and
`cached_friedman_results.sig_by_location` /
`suppression_by_location`:

| # | Feature | Source | Description |
|---|---------|--------|-------------|
| 1 | `contra_vis_activation` | group_mean_frs | Mean (visual FR - baseline FR) at contralateral locations |
| 2 | `ipsi_vis_activation` | group_mean_frs | Mean (visual FR - baseline FR) at ipsilateral locations |
| 3 | `contra_sac_activation` | group_mean_frs | Mean (saccade FR - baseline FR) at contralateral locations |
| 4 | `ipsi_sac_activation` | group_mean_frs | Mean (saccade FR - baseline FR) at ipsilateral locations |
| 5 | `contra_delay_activation` | group_mean_frs | Mean (delay FR - baseline FR) at contralateral locations |
| 6 | `ipsi_delay_activation` | group_mean_frs | Mean (delay FR - baseline FR) at ipsilateral locations |
| 7 | `spatial_selectivity_vis` | derived | (contra_vis - ipsi_vis) / (|contra_vis| + |ipsi_vis| + eps) |
| 8 | `spatial_selectivity_sac` | derived | (contra_sac - ipsi_sac) / (|contra_sac| + |ipsi_sac| + eps) |
| 9 | `n_locations_activated` | sig_by_location | Count of locations with significant activation (any epoch) |
| 10 | `n_locations_suppressed` | suppression_by_location | Count of locations with significant suppression |
| 11 | `vis_response_magnitude` | group_mean_frs | Max |FR change| in visual epoch across all locations |
| 12 | `sac_response_magnitude` | group_mean_frs | Max |FR change| in saccade epoch across all locations |
| 13 | `delay_response_magnitude` | group_mean_frs | Max |FR change| in delay epoch across all locations |
| 14 | `vis_to_sac_ratio` | derived | vis_magnitude / (vis_magnitude + sac_magnitude + eps) |
| 15 | `delay_index` | derived | delay_magnitude / (max(vis, sac) + eps) |
| 16 | `has_suppression` | suppression_by_location | Any suppression anywhere (0/1) |
| 17 | `has_mixed_polarity` | both | Both activation AND suppression present (0/1) |

**Key implementation details:**

1. Identify contra/ipsi locations using `location_angles` and
   `scSide` (from screen_sc_neurons). Need scSide - get it from
   the cached results or re-derive from target angles.

2. For gSac_jph group (last group if present): treat as
   contralateral (it's always placed at RF center).

3. FR changes are computed as (epoch_mean - baseline_mean) from
   `group_mean_frs{i_group}(i_cluster, :)` where columns are
   [baseline, visual, delay, saccade].

4. Only compute features for neurons passing quality gates.
   Return NaN rows for excluded neurons.

5. Normalize features using z-score (`zscore`) across the
   included population ONLY.

**Contra/Ipsi determination:**
```matlab
% From screen_sc_neurons, scSide is determined.
% We need the same logic here.
% Right SC: contralateral = left VF (90 < theta < 270)
% Left SC: contralateral = right VF (theta < 90 | theta > 270)
% Use location_angles from cached_friedman_results
```

---

## 3. Phase 2: Clustering Pipeline

### File: `code/cluster_sc_neurons.m`

**Purpose:** Determine optimal clusters and assign neurons.

**Signature:**
```matlab
function clustering = cluster_sc_neurons( ...
    feature_matrix, feature_names, cluster_ids)
% Returns: clustering struct with fields:
%   .assignments     [N x 1] cluster index (1..k)
%   .k               optimal number of clusters
%   .centroids       [k x F] cluster centroids
%   .silhouette_avg  mean silhouette score
%   .silhouette_per_neuron [N x 1]
%   .evaluation      struct with BIC, silhouette for each k
%   .pca_scores      [N x nPC] PCA-transformed features
%   .pca_explained   [nPC x 1] variance explained
%   .pca_coeff       [F x nPC] PCA loadings
%   .labels          {k x 1} cluster labels (initially
%                    'cluster_1', ..., assigned later)
```

**Algorithm:**

1. **PCA** (optional, for visualization and correlated-feature removal):
   ```matlab
   [coeff, scores, ~, ~, explained] = pca(feature_matrix);
   n_pc = find(cumsum(explained) >= 90, 1);
   pca_features = scores(:, 1:n_pc);
   ```
   Use PCA-transformed features for clustering ONLY if many
   features are highly correlated (r > 0.8). Otherwise cluster
   on z-scored raw features.

2. **Determine optimal k** using `evalclusters`:
   ```matlab
   % Test k = 2 to min(10, floor(N/15))
   % floor(N/15) ensures at least 15 neurons per cluster
   max_k = min(10, floor(N / 15));

   eva_sil = evalclusters(features, 'kmeans', ...
       'silhouette', 'KList', 2:max_k);
   eva_gap = evalclusters(features, 'kmeans', ...
       'gap', 'KList', 2:max_k);

   % Also compute BIC for GMM
   bic_vals = nan(1, max_k);
   for k = 2:max_k
       gm = fitgmdist(features, k, ...
           'RegularizationValue', 0.01, ...
           'Replicates', 5);
       bic_vals(k) = gm.BIC;
   end

   % Choose k: prefer silhouette optimal, cross-check
   % with gap and BIC
   optimal_k = eva_sil.OptimalK;
   ```

3. **Final clustering** with k-means:
   ```matlab
   rng(42);  % Fixed seed for reproducibility
   [idx, centroids] = kmeans(features, optimal_k, ...
       'Replicates', 20, ...
       'MaxIter', 500, ...
       'Distance', 'sqeuclidean');
   ```

4. **Also run hierarchical clustering** for comparison:
   ```matlab
   Z = linkage(features, 'ward');
   idx_hier = cluster(Z, 'maxclust', optimal_k);
   % Compare to k-means via adjusted Rand index
   ```

5. **Compute silhouette scores:**
   ```matlab
   [sil_vals, ~] = silhouette(features, idx);
   mean_sil = mean(sil_vals);
   ```

6. **Cluster stability** (bootstrap):
   ```matlab
   n_boot = 100;
   stability = nan(N, n_boot);
   for b = 1:n_boot
       boot_idx = randsample(N, N, true);
       [boot_assignments, ~] = kmeans( ...
           features(boot_idx, :), optimal_k, ...
           'Start', centroids, 'MaxIter', 100);
       % Map boot clusters to original via Hungarian
       stability(:, b) = mapped_assignments;
   end
   % Proportion of bootstrap runs agreeing with
   % original assignment
   cluster_stability = mean(stability == idx, 2);
   ```

---

## 4. Phase 3: Visualization & Interpretation

### File: `code/visualize_sc_clusters.m`

**Purpose:** Generate diagnostic figures for cluster interpretation.

**Signature:**
```matlab
function visualize_sc_clusters(features, clustering, ...
    v1_classes, project_root)
```

**Figures to generate:**

1. **PCA scatter** (PC1 vs PC2, colored by cluster):
   - Also show v1 classification with different markers
   - Saved as `figures/sc_clustering_pca.pdf`

2. **Feature radar plots** (one per cluster):
   - Mean feature values per cluster on a radar/spider chart
   - Identifies dominant characteristics
   - Saved as `figures/sc_clustering_radar.pdf`

3. **Silhouette plot:**
   - Standard silhouette diagram grouped by cluster
   - Saved as `figures/sc_clustering_silhouette.pdf`

4. **Cluster evaluation plot:**
   - Silhouette vs k, gap vs k, BIC vs k
   - Saved as `figures/sc_clustering_evaluation.pdf`

5. **Confusion matrix: v1 vs v2:**
   - Crosstab of old classification vs new clusters
   - Shows how classic/interneuron map to new categories
   - Saved as `figures/sc_clustering_v1_vs_v2.pdf`

6. **Dendrogram** (from hierarchical clustering):
   - Color-coded by final cluster assignment
   - Saved as `figures/sc_clustering_dendrogram.pdf`

---

## 5. Phase 4: Validation

### File: `code/validate_sc_classification.m`

**Purpose:** Assess classification quality.

**Metrics computed:**

1. **Silhouette score**: Target > 0.4
2. **Cluster stability** (bootstrap): Target > 80% assignment
   agreement
3. **Cross-session stability**: Leave-one-session-out, train on
   N-1 sessions, classify held-out session using nearest
   centroid. Report agreement rate.
4. **v1 vs v2 concordance**: What fraction of v1 "classic"
   neurons land in the v2 "classic" cluster? (Expect high overlap)

**Output:** Prints summary table and saves
`figures/sc_classification_validation.pdf`

**Expert validation** is deferred to runtime: after running the
pipeline, the user reviews exemplar neurons per cluster using
existing neuron summary PDFs.

---

## 6. Phase 5: Pipeline Integration

### 6.1 Modify `screen_sc_neurons.m`

**Changes (minimal, additive):**

After the existing classification block (line 636), add:

```matlab
%% Data-Driven Classification (v2)
% Extract features and cluster neurons that passed
% quality gates
features = extract_sc_response_features(session_data);
if sum(features.passes_quality) >= 30
    clustering = cluster_sc_neurons( ...
        features.feature_matrix(features.passes_quality, :), ...
        features.feature_names, ...
        features.cluster_ids(features.passes_quality));

    % Map cluster assignments back to full neuron array
    neuron_class_v2 = repmat({'excluded'}, nClusters, 1);
    quality_idx = find(features.passes_quality);
    for i = 1:length(quality_idx)
        neuron_class_v2{quality_idx(i)} = ...
            clustering.labels{clustering.assignments(i)};
    end
else
    % Too few neurons for reliable clustering
    neuron_class_v2 = neuron_class;  % Fall back to v1
    clustering = [];
end
```

**New output:** Add `neuron_class_v2` as 6th output argument.

**Updated signature:**
```matlab
function [selected_neurons, sig_epoch_comparison, scSide, ...
    neuron_class, session_data, neuron_class_v2] = ...
    screen_sc_neurons(session_data, project_root, ...
    condition_defs, session_name)
```

### 6.2 Modify `apply_neuron_screening.m`

**Changes:**

In the SC-specific classification block (~line 384-417):

```matlab
% After existing code that sets screening_info(i).neuron_class:
% Add v2 classification
if exist('sc_neuron_class_v2', 'var')
    for i = 1:nNeurons
        if i <= length(sc_neuron_class_v2)
            screening_info(i).neuron_class_v2 = ...
                sc_neuron_class_v2{i};
        else
            screening_info(i).neuron_class_v2 = 'excluded';
        end
    end
end
```

Also add `neuron_class_v2` to the screening_info initialization
(~line 154):
```matlab
screening_info(i).neuron_class_v2 = 'excluded';
```

### 6.3 Modify `aggregate_analysis_results.m`

**Changes:**

In the neuron_type extraction block (~line 512-546), add
parallel extraction of v2:

```matlab
neuron_type_v2_col = repmat({''}, n_session_neurons, 1);
for i_neuron = 1:n_session_neurons
    if isfield(screening_info(i_neuron), 'neuron_class_v2')
        neuron_type_v2_col{i_neuron} = ...
            screening_info(i_neuron).neuron_class_v2;
    end
end
```

Add `neuron_type_v2_col` and `'neuron_type_v2'` to the session
table construction (~line 606-610).

### 6.4 Downstream plotting (deferred)

`plot_aggregated_window_roc.m` currently has hardcoded
`'classic'` and `'interneuron'` checks (lines 113-114). This
will need updating once v2 is validated and promoted, but
NOT during initial implementation. The v1 plotting continues
to work unchanged.

---

## 7. File Summary

### New Files

| File | Phase | Purpose |
|------|-------|---------|
| `code/extract_sc_response_features.m` | 1 | Feature extraction from cached_friedman_results |
| `code/cluster_sc_neurons.m` | 2 | Clustering pipeline (optimal k, k-means, stability) |
| `code/visualize_sc_clusters.m` | 3 | Cluster diagnostic figures |
| `code/validate_sc_classification.m` | 4 | Validation metrics |

### Modified Files

| File | Phase | Changes |
|------|-------|---------|
| `code/screen_sc_neurons.m` | 5 | Add v2 classification block, new output arg |
| `code/apply_neuron_screening.m` | 5 | Propagate neuron_class_v2 |
| `code/aggregate_analysis_results.m` | 5 | Add neuron_type_v2 column |
| `docs/design/analysis_rationale.md` | 0 | Already updated (Section 3.2) |

### NOT Modified (yet)

| File | Reason |
|------|--------|
| `code/plot_aggregated_window_roc.m` | Deferred until v2 validated |
| `code/define_task_conditions.m` | No threshold changes needed |
| `code/run_4factors_analysis.m` | No changes needed (passes through) |

---

## 8. Implementation Order

```
Step 1: extract_sc_response_features.m
        - Can be tested standalone on any session_data
        - Verify feature distributions look reasonable

Step 2: cluster_sc_neurons.m
        - Requires aggregated features across all sessions
        - Need to run feature extraction on ALL sessions first
        - Build a driver script to aggregate features and cluster

Step 3: visualize_sc_clusters.m
        - Generates diagnostic figures for expert review
        - User reviews and provides cluster label feedback

Step 4: validate_sc_classification.m
        - Computes quality metrics
        - User reviews exemplar neurons per cluster

Step 5: Integration (screen_sc_neurons.m changes)
        - Only after Steps 3-4 confirm quality

Step 6: Propagation (apply_neuron_screening.m,
        aggregate_analysis_results.m)
        - Only after Step 5
```

**Critical dependency:** Steps 1-4 can be developed and tested
independently of the main pipeline. Step 5-6 modifies the
pipeline and should only happen after validation.

For Steps 1-4, we need a **driver script** that:
1. Loads all SC sessions
2. Runs feature extraction on each
3. Aggregates feature matrices
4. Runs clustering on the aggregate
5. Generates visualizations

### File: `code/run_sc_clustering_analysis.m` [NEW]

This is a standalone script (not called by the main pipeline)
that runs the full clustering analysis for development and
validation.

---

## 9. Coding Conventions Checklist

- [ ] Max 75 chars/line with `...` continuation
- [ ] Comments above code, not to the right
- [ ] `mySubPlot([a,b,c])` not `subplot(a,b,c)`
- [ ] `barStairsFill.m` for any PSTHs
- [ ] `richColors.m` for color palettes
- [ ] `'Interpreter','none'` for titles with underscores
- [ ] Path via `findOneDrive` pattern
- [ ] `%%` section markers for organization
- [ ] Function headers with INPUTS/OUTPUTS documentation

---

## 10. Testing Strategy

1. **Unit test feature extraction**: Run on a single session,
   verify feature values match manual calculation from
   cached_friedman_results.

2. **Sanity check clustering**: With k=2, verify it roughly
   separates classic from interneuron (v1 categories).

3. **Visual inspection**: Generate PCA scatter colored by v1
   classification. Expect classic neurons to cluster together.

4. **Smoke test pipeline**: After integration, run full pipeline
   on one session and verify no errors.

5. **Regression test**: Compare v1 neuron_class output before
   and after changes - must be identical.
