# SC Classification v2 - Feature Redesign Plan

## Background

The initial v2 implementation used 17 features that collapsed 4 target
locations into contra/ipsi averages. The user rejected this design:
features should use all 4 locations independently. The first clustering
run also produced k=2 (13 vs 416 neurons), an unhelpfully imbalanced
result, motivating a minimum cluster size constraint.

## Design Decision: Rank-Ordered Per-Location Features

**Problem:** Target angles vary across sessions (45 vs 67 vs 90 degrees
for "location 1"). Using raw location indices would conflate neuron
type with RF position.

**Solution:** For each neuron, sort the 4 gSac_4factors locations by
visual response magnitude (baseline-subtracted). Features become
best/2nd/3rd/worst, making them RF-position invariant. Two neurons
with identical response profiles at different RF positions get
identical feature vectors.

## New Feature Set (19 features)

### Per-Location FR Changes (10 features)
Locations ranked by visual response (vis - baseline), descending:

| # | Feature | Description |
|---|---------|-------------|
| 1 | `best_vis` | FR change at best visual location |
| 2 | `2nd_vis` | FR change at 2nd-best visual location |
| 3 | `3rd_vis` | FR change at 3rd-best visual location |
| 4 | `worst_vis` | FR change at worst visual location |
| 5 | `best_sac` | Saccade FR change at best-visual location |
| 6 | `2nd_sac` | Saccade FR change at 2nd-best-visual loc |
| 7 | `3rd_sac` | Saccade FR change at 3rd-best-visual loc |
| 8 | `worst_sac` | Saccade FR change at worst-visual loc |
| 9 | `best_delay` | Delay FR change at best-visual location |
| 10 | `worst_delay` | Delay FR change at worst-visual location |

Note: sac and delay use the same rank order as visual. This preserves
the spatial relationship (does the best visual location also have a
saccade response?).

### Significance Counts (4 features)

| # | Feature | Description |
|---|---------|-------------|
| 11 | `n_loc_sig_vis` | Locations with significant visual activation |
| 12 | `n_loc_sig_sac` | Locations with significant saccade activation |
| 13 | `n_loc_sig_delay` | Locations with significant delay activation |
| 14 | `n_loc_suppressed` | Locations with any suppression |

### Temporal Ratios (2 features)

| # | Feature | Description |
|---|---------|-------------|
| 15 | `vis_to_sac_ratio` | best_vis / (best_vis + best_sac + eps) |
| 16 | `delay_index` | best_delay / max(best_vis, best_sac) + eps |

### Polarity Features (2 features)

| # | Feature | Description |
|---|---------|-------------|
| 17 | `has_suppression` | Any suppression at any location/epoch (0/1) |
| 18 | `has_mixed_polarity` | Both activation AND suppression present (0/1) |

### Tuning Features (1 feature)

| # | Feature | Description |
|---|---------|-------------|
| 19 | `tuning_breadth` | best_vis / mean(all 4 vis changes + eps) |

`tuning_breadth` captures how sharply tuned the neuron is. A value
near 1 means all locations respond equally (omnidirectional). A high
value means only one location responds strongly (narrow RF).

## gSac_jph Handling

The 5th location group (gSac_jph trials) is excluded from feature
extraction. It appears in only 4/7 sessions, is always contralateral,
and represents a different task condition. Only the 4 gSac_4factors
locations are used.

## Files to Modify

### 1. `code/extract_sc_response_features.m` (REWRITE)
- Remove contra/ipsi logic entirely
- Filter to first 4 location groups (gSac_4factors only)
- Compute FR change from baseline at each location x 3 epochs
- Sort locations by visual FR change (descending) per neuron
- Extract 19 features using rank-ordered data
- Quality gate logic unchanged

### 2. `code/cluster_sc_neurons.m` (ADD MIN CLUSTER SIZE)
- After k-means, check if any cluster has < `min_per_cluster` neurons
- If so, merge small cluster with nearest centroid neighbor
- Re-compute silhouette after merge
- Report merges in console output

### 3. `code/run_sc_clustering_analysis.m` (UPDATE FEATURE NAMES)
- Update hardcoded `feat_names` cell array to match new 19 features
- `nFeatures` updates automatically

### 4. `code/visualize_sc_clusters.m` (NO CHANGES NEEDED)
- Feature bar plots use `features_agg.feature_names` dynamically
- Will auto-adapt to new feature set

## Implementation Order

1. Rewrite `extract_sc_response_features.m`
2. Add min cluster size to `cluster_sc_neurons.m`
3. Update `run_sc_clustering_analysis.m` feature names
4. Run pipeline end-to-end
5. Review results

## Verification Criteria

- All 7 SC sessions load and extract features successfully
- Feature matrix has 19 columns, no all-NaN columns
- Quality-passing neuron count matches previous run (~429)
- k > 2 (improved over previous run)
- No cluster has fewer than 15 neurons
- Silhouette scores printed and saved
- All 6 diagnostic figures generated
- CSV export and exemplar neuron printout work correctly
