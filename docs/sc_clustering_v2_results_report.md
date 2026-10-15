# SC Classification v2 - Results Report

**Date:** 2026-02-06
**Author:** Claude Code

## Summary

The SC neuron classification v2 pipeline was redesigned and executed.
The key change from the initial implementation was replacing contra/ipsi-
averaged features with rank-ordered per-location features (19 total),
where all 4 target locations contribute independently. The clustering
identified **k=4** neuron groups from 411 quality-passing neurons
across 7 SC sessions.

## What Changed from v1 Implementation

### Feature Extraction (`extract_sc_response_features.m`)
- **Old:** 17 features, collapsed 4 locations into contra/ipsi averages
- **New:** 19 features, rank-ordered by visual response strength
  (best/2nd/3rd/worst), making features RF-position invariant
- gSac_jph (5th location group) excluded from features
- See `docs/sc_clustering_v2_feature_redesign.md` for full spec

### Clustering (`cluster_sc_neurons.m`)
- Added minimum cluster size enforcement (merges small clusters)
- Default min_per_cluster lowered from 15 to 10 (12-neuron cluster
  was genuinely distinct; merging it degraded overall quality)

### Driver (`run_sc_clustering_analysis.m`)
- Feature names updated to match new extraction

## Results

### Pipeline Execution
- 7 SC sessions processed
- 694 total neurons, 411 pass quality gates (59.2%)
- 19 features extracted per neuron
- 11 PCA components capture 92.1% variance
- k selection: Silhouette optimal=4, Gap optimal=10, BIC optimal=4
- Selected k=4 (silhouette + BIC consensus)
- No clusters merged (all >= 10 neurons)

### Cluster Summary

| Cluster | N | % | Silhouette | Bootstrap | v1 Classic | v1 Interneuron | v1 Excluded |
|---------|---|---|-----------|-----------|-----------|---------------|-------------|
| 1 | 113 | 27.5% | -0.103 | 98.8% | 26 | 57 | 30 |
| 2 | 21 | 5.1% | -0.142 | 85.9% | 10 | 11 | 0 |
| 3 | 265 | 64.5% | 0.676 | 99.9% | 89 | 52 | 124 |
| 4 | 12 | 2.9% | 0.304 | 99.5% | 2 | 8 | 2 |
| **Mean** | | | **0.409** | **98.9%** | | | |

### Quality Metrics
- Mean silhouette: 0.409 (PASS, threshold > 0.4)
- Bootstrap stability: 98.9% (PASS, threshold > 80%)
- Leave-one-session-out agreement: 95.5% (PASS, threshold > 70%)

### Cluster Interpretation (Preliminary)

**Cluster 3 (265 neurons, 64.5%):** The dominant group. Well-separated
(silhouette 0.676). Enriched for v1 "classic" neurons (89 of 127 total
classic neurons = 70%). Also contains 124 v1-excluded neurons that
pass v2 quality gates, suggesting v1 may have been too conservative
in excluding modulated neurons lacking spatial selectivity. These
neurons likely have modest but consistent responses.

**Cluster 1 (113 neurons, 27.5%):** Interneuron-enriched (57 of 128
total interneurons = 44.5%). Negative silhouette (-0.103) indicates
overlap with other clusters in feature space. These likely represent
neurons with suppressed or mixed-polarity responses.

**Cluster 4 (12 neurons, 2.9%):** Small but genuinely distinct
(silhouette 0.304). Strongly interneuron-enriched (8 of 12 from v1
interneurons). May represent a specific interneuron subtype with
distinctive response properties.

**Cluster 2 (21 neurons, 5.1%):** Mixed classic/interneuron composition.
Negative silhouette (-0.142) and lower bootstrap stability (85.9%).
This is the least well-defined cluster — may represent a transition
zone between types.

## Quality Concerns

1. **Negative silhouettes in clusters 1 and 2:** These neurons are on
   average closer to a neighboring cluster than their own. This could
   indicate (a) continuous rather than discrete neuron types, or
   (b) the feature set doesn't fully capture the relevant distinctions.
   This is common in neural data — biological categories often have
   fuzzy boundaries.

2. **Cluster 2 is small (21 neurons):** Borderline for statistical
   reliability. Its distinctive features should be interpreted
   cautiously.

3. **Cluster 3 dominance (64.5%):** A single cluster containing nearly
   two-thirds of neurons suggests the primary split is "responsive vs
   differently-responsive" rather than multiple well-separated types.

## Output Files

| File | Description |
|------|-------------|
| `figures/sc_clustering_pca.pdf` | PCA scatter: v2 clusters and v1 classes |
| `figures/sc_clustering_features.pdf` | Mean feature values per cluster with SEM |
| `figures/sc_clustering_silhouette.pdf` | Silhouette plot |
| `figures/sc_clustering_evaluation.pdf` | k selection curves (silhouette, gap, BIC) |
| `figures/sc_clustering_v1_vs_v2.pdf` | v1 vs v2 confusion matrix |
| `figures/sc_clustering_dendrogram.pdf` | Hierarchical clustering dendrogram |
| `figures/sc_clustering_neuron_assignments.csv` | Per-neuron cluster assignments |
| `data/sc_clustering_results.mat` | Full results (clustering, validation, features) |

## Exemplar Neurons (Closest to Centroid)

### Cluster 1 (113 neurons)
1. Feynman_08_13_2025_SC cluster 50 (v1: excluded)
2. Feynman_08_13_2025_SC cluster 72 (v1: excluded)
3. Newton_08_13_2025_SC cluster 32 (v1: excluded)
4. Feynman_08_12_2025_SC cluster 0 (v1: classic)
5. Feynman_08_13_2025_SC cluster 64 (v1: excluded)

### Cluster 2 (21 neurons)
1. Newton_08_13_2025_SC cluster 84 (v1: classic)
2. Newton_08_15_2025_SC cluster 76 (v1: classic)
3. Newton_08_13_2025_SC cluster 72 (v1: classic)
4. Newton_08_15_2025_SC cluster 67 (v1: interneuron)
5. Newton_08_15_2025_SC cluster 85 (v1: classic)

### Cluster 3 (265 neurons)
1. Newton_08_13_2025_SC cluster 86 (v1: excluded)
2. Feynman_08_15_2025_SC cluster 63 (v1: classic)
3. Newton_08_15_2025_SC cluster 3 (v1: classic)
4. Feynman_08_12_2025_SC cluster 6 (v1: excluded)
5. Feynman_08_12_2025_SC cluster 18 (v1: classic)

### Cluster 4 (12 neurons)
1. Newton_08_14_2025_SC cluster 48 (v1: excluded)
2. Newton_08_14_2025_SC cluster 76 (v1: excluded)
3. Newton_08_15_2025_SC cluster 84 (v1: classic)
4. Newton_08_13_2025_SC cluster 71 (v1: interneuron)
5. Newton_08_14_2025_SC cluster 62 (v1: interneuron)

## Next Steps

1. **Review diagnostic PDFs** for the exemplar neurons listed above.
   Open the corresponding `*_neuron_diagnostics.pdf` files in the
   session figure directories.

2. **Review the feature distribution PDF** (`sc_clustering_features.pdf`)
   to understand what distinguishes each cluster.

3. **Open `sc_clustering_neuron_assignments.csv`** to browse all
   neuron assignments. Columns include session_id, cluster_id,
   v2_cluster, distance_to_centroid, v1_class, silhouette, and
   boot_stability.

4. **Assign meaningful labels** to the 4 clusters based on your
   neuroscientific expertise. The current labels are `cluster_1`
   through `cluster_4`. Update them in the .mat file:
   ```matlab
   load('data/sc_clustering_results.mat');
   clustering.labels = {'label1', 'label2', ...
       'label3', 'label4'};
   save('data/sc_clustering_results.mat', ...
       'clustering', '-append');
   ```

5. **Consider whether k=4 is appropriate.** If the clusters don't
   make neuroscientific sense, you can re-run with `force_k` set
   to a different value, or adjust min_per_cluster.

## Code Files Modified

| File | Change |
|------|--------|
| `code/extract_sc_response_features.m` | Rewritten: rank-ordered per-location features |
| `code/cluster_sc_neurons.m` | Added min cluster size enforcement, default=10 |
| `code/run_sc_clustering_analysis.m` | Updated feature names to match new extraction |
| `code/visualize_sc_clusters.m` | No changes (auto-adapts to feature names) |
| `code/validate_sc_classification.m` | No changes |
| `code/screen_sc_neurons.m` | No changes (from iteration 1) |
| `code/apply_neuron_screening.m` | No changes (from iteration 1) |
| `code/aggregate_analysis_results.m` | No changes (from iteration 1) |
