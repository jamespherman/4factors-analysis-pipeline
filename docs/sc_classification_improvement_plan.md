# SC Neuron Classification Improvement Plan

**Created:** 2026-02-04
**Status:** Planning document - awaiting review and implementation

---

## Executive Summary

This plan addresses the need to improve SC neuron classification beyond the current ad-hoc "classic" vs "interneuron" distinction. The goal is to leverage the unique multi-location response data from the 4factors task to:

1. Characterize SC neurons using their full spatial response profile (4 locations)
2. Use data-driven clustering to identify natural groupings
3. Map clusters to biologically meaningful categories
4. Create a reproducible, algorithmic classification that matches expert intuition

**Key insight:** Traditional SC classification used single-location RF mapping. The 4factors task provides responses at 4 locations (2 contra, 2 ipsi), enabling characterization of neurons that would have been "skipped over" in classical studies due to lack of clear spatial tuning.

---

## 1. Current State Analysis

### 1.1 Current Classification Logic

From `screen_sc_neurons.m`, the current approach:

**Classic SC neuron** (all must be true):
- Significant activation vs baseline at ≥1 contralateral location (any epoch)
- Mean FR > 3.0 sp/s
- NOT significantly activated at opposing locations (180° apart) during Visual/Saccade epochs
- Not sparse (≤80% empty bins)

**Putative interneuron:**
- Mean FR > 3.0 sp/s
- Task-modulated somewhere
- Fails "classic" criteria (has opposing activation OR contralateral suppression)

### 1.2 Current Results

| Category | Count | % of Total |
|----------|-------|------------|
| Classic | 138 | 19.9% |
| Interneuron | 131 | 18.9% |
| Excluded | 425 | 61.2% |
| **Total SC** | **694** | 100% |

### 1.3 Identified Problems

1. **Sensitivity to weak ipsilateral responses**: Neurons with strong contralateral responses but statistically significant (though small magnitude) ipsilateral responses get classified as "interneuron" due to the "opposing activation" criterion.

2. **Heterogeneous "interneuron" category**: This category conflates:
   - Broadly activated neurons (respond to multiple locations)
   - Suppressed neurons (contralateral suppression)
   - Neurons lacking contralateral preference

3. **Binary spatial criterion**: The current approach only checks whether opposing locations are BOTH activated, not the relative magnitude or pattern.

4. **No consideration of epoch-specific patterns**: Traditional V/VM/M/D classification (Visual/Visuomotor/Motor/Delay) is not captured.

---

## 2. Literature Context

### 2.1 Traditional SC Classification (Single-Location)

The classic approach ([Phongphanphanee et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37202508/)) classifies neurons into:

| Class | Visual Epoch | Delay Epoch | Motor Epoch |
|-------|--------------|-------------|-------------|
| Visual (V) | ++ | - | - |
| Visuomotor (VM) | + | +/- | + |
| Motor (M) | - | - | ++ |
| Delay (D) | +/- | ++ | +/- |

This classification assumes a single contralateral target location within the RF.

### 2.2 Inhibitory Neurons in SC

Research on GABAergic neurons in SC ([Whyland et al., 2023](https://www.frontiersin.org/articles/10.3389/fncir.2018.00035); [Lee & bharioke, 2021](https://www.nature.com/articles/s42003-021-02248-1)) reveals:

- **High density**: SC has among the highest densities of inhibitory neurons in the brain
- **Morphological types**: Intralaminar interneurons (local/horizontal), interlaminar interneurons, commissural neurons
- **Functional roles**: Local inhibition, inter-hemispheric inhibition, gain control
- **Spatial properties**: Some GABAergic neurons project to distal SC regions (same or opposite SC)

**Key insight for classification:** True inhibitory interneurons may show:
- Suppression in response to contralateral stimuli (if they receive input from excitatory neurons and provide feedback inhibition)
- Activation that is spatially broad (if they pool across locations)
- Activation triggered by ipsilateral stimuli (if they receive commissural input)

### 2.3 Data-Driven Classification Approaches

Modern approaches use unsupervised clustering ([Goetz et al., 2022](https://elifesciences.org/articles/82367); [Automatic neuron classification](https://pmc.ncbi.nlm.nih.gov/articles/PMC4417416/)):

- **Feature spaces**: Morphology, physiology, molecular markers
- **Algorithms**: Hierarchical clustering, k-means, Gaussian mixture models
- **Evaluation**: Bayesian Information Criterion (BIC) for optimal cluster number
- **Handling overlap**: Probabilistic assignment, soft clustering

---

## 3. Proposed Approach

### 3.1 Philosophy

Rather than defining categories a priori, we will:

1. **Extract a rich feature set** from each neuron's multi-location responses
2. **Apply unsupervised clustering** to find natural groupings
3. **Interpret clusters** in terms of biologically meaningful categories
4. **Validate** against expert annotation and functional predictions

### 3.2 Feature Extraction

For each neuron, extract the following features from the cached Friedman test results:

#### A. Spatial Response Profile (per epoch: Visual, Delay, Saccade)

| Feature | Description | Range |
|---------|-------------|-------|
| `contra_activation_strength` | Mean FR increase at contralateral locations | 0 to max |
| `ipsi_activation_strength` | Mean FR increase at ipsilateral locations | 0 to max |
| `contra_suppression_strength` | Mean FR decrease at contralateral locations | 0 to max |
| `ipsi_suppression_strength` | Mean FR decrease at ipsilateral locations | 0 to max |
| `spatial_selectivity_index` | (Contra - Ipsi) / (Contra + Ipsi) | -1 to +1 |
| `n_locations_activated` | Count of locations with significant activation | 0 to 4 |
| `n_locations_suppressed` | Count of locations with significant suppression | 0 to 4 |

#### B. Temporal Response Profile

| Feature | Description |
|---------|-------------|
| `visual_response` | Max FR change in visual epoch (any location) |
| `delay_activity` | Max FR change in delay epoch |
| `saccade_response` | Max FR change in saccade epoch |
| `vis_to_sac_ratio` | Visual response / Saccade response |
| `delay_index` | Delay activity / max(Visual, Saccade) |

#### C. Response Polarity

| Feature | Description |
|---------|-------------|
| `is_excited_anywhere` | Any location shows activation |
| `is_suppressed_anywhere` | Any location shows suppression |
| `mixed_polarity` | Both excitation and suppression present |

### 3.3 Clustering Pipeline

```
Step 1: Feature extraction
    └── Extract features for all neurons passing basic criteria
        (FR > threshold, not sparse)

Step 2: Feature normalization
    └── Z-score normalize each feature across neurons

Step 3: Dimensionality reduction (optional)
    └── PCA to reduce to top N components capturing 90% variance
    └── Enables visualization and removes correlated features

Step 4: Determine optimal cluster number
    └── Apply hierarchical clustering with Ward's linkage
    └── Compute silhouette scores for k = 2 to 10
    └── Use elbow method and BIC to select k

Step 5: Cluster assignment
    └── Apply k-means or hierarchical clustering with chosen k
    └── Compute cluster centroids and within-cluster variance

Step 6: Cluster interpretation
    └── Analyze feature distributions within each cluster
    └── Map to biologically meaningful categories
    └── Generate interpretable cluster labels
```

### 3.4 Proposed Category Framework

Based on literature and the user's observations, we anticipate the following categories may emerge:

| Category | Expected Features | Biological Interpretation |
|----------|-------------------|--------------------------|
| **Classic-Visual** | Strong contra activation (visual epoch), weak ipsi, no suppression | Traditional SC visual neurons |
| **Classic-Motor** | Strong contra activation (saccade epoch), weaker visual | Traditional SC motor neurons |
| **Classic-Visuomotor** | Strong contra activation (both visual and saccade epochs) | Traditional SC visuomotor |
| **Delay-Active** | Prominent delay epoch activity, spatially selective | Memory-related activity |
| **Broad-Activated** | Activation at multiple locations (including ipsi), not spatially selective | Potential wide-field neurons or pooling interneurons |
| **Suppressed** | Contralateral suppression, may have ipsi activation | Potential inhibitory interneurons |
| **Mixed-Polarity** | Both activation and suppression at different locations | Complex spatial integration |

**Important:** These are hypothesized categories. The actual clusters may differ, and the final labels should reflect the data.

---

## 4. Implementation Plan

### Phase 1: Data Extraction and Exploration (2-3 hours)

**Step 1.1: Create feature extraction function**

```matlab
function features = extract_sc_response_features(session_data)
% EXTRACT_SC_RESPONSE_FEATURES Extract multi-location response features.
%
% Uses cached_friedman_results to extract spatial and temporal response
% features for each neuron.
%
% OUTPUT: features table with one row per neuron
```

**Step 1.2: Aggregate features across sessions**

- Load all SC sessions
- Extract features for neurons passing basic criteria
- Combine into single feature matrix

**Step 1.3: Exploratory visualization**

- Feature histograms
- Pairwise scatter plots (e.g., contra vs ipsi activation)
- PCA scatter plot colored by current classification

### Phase 2: Clustering Analysis (2-3 hours)

**Step 2.1: Determine feature set**

- Compute feature correlations
- Remove highly correlated features (r > 0.8)
- Apply PCA if needed

**Step 2.2: Cluster number selection**

- Hierarchical clustering dendrogram
- Silhouette analysis for k = 2 to 10
- BIC/AIC for mixture models

**Step 2.3: Final clustering**

- Apply chosen algorithm with optimal k
- Compute cluster stability (bootstrap resampling)
- Assign each neuron to a cluster

### Phase 3: Cluster Interpretation (1-2 hours)

**Step 3.1: Analyze cluster centroids**

- Generate radar plots of mean features per cluster
- Compute statistical comparisons between clusters

**Step 3.2: Map to categories**

- Review cluster characteristics
- Assign biologically meaningful labels
- Document rationale for each label

**Step 3.3: Generate example neurons**

- Select 3-5 representative neurons per cluster
- Review their individual summary PDFs
- Validate cluster assignment makes sense

### Phase 4: Validation (2-3 hours)

**Step 4.1: Expert annotation**

- Select ~30 neurons (5 per anticipated cluster)
- User classifies based on visual inspection of PDFs
- Compare to algorithmic classification
- Compute agreement rate (target: >80%)

**Step 4.2: Cross-validation**

- Leave-one-session-out validation
- Train clustering on N-1 sessions
- Classify held-out session
- Assess cluster stability across sessions

**Step 4.3: Functional predictions**

- Test whether cluster membership predicts 4-factors selectivity
- Compare to current classification's predictive power

### Phase 5: Integration (1-2 hours)

**Step 5.1: Update screening function**

- Modify `screen_sc_neurons.m` to use new classification
- Maintain backward compatibility (add new fields, don't remove old)

**Step 5.2: Update aggregation**

- Add new `neuron_class_detailed` field
- Preserve original `neuron_class` for comparison

**Step 5.3: Documentation**

- Update `analysis_rationale.md` with new classification
- Document cluster characteristics and rationale
- Add to PROGRESS_AND_TODO.md

---

## 5. Specific Questions to Address

During implementation, we should investigate:

### Q1: Magnitude vs Statistical Significance

**Issue:** Current approach uses only significance (p < alpha). A neuron with FR change from 5 to 5.5 sp/s could be "significant" with enough trials, while one with FR change from 5 to 15 sp/s might not be if variable.

**Question:** Should classification use:
- Only significance (current)
- Only magnitude (e.g., >50% FR increase)
- Both (significant AND >X% change)
- Effect size (Cohen's d or similar)

**Recommendation:** Use a combined criterion: significant (p < 0.05/3) AND magnitude > 20% change from baseline.

### Q2: Epoch Weighting

**Issue:** Current approach treats Visual, Delay, and Saccade epochs equally for spatial selectivity but only considers Visual/Saccade for "opposing activation."

**Question:** Should different epochs be weighted differently? Should spatial selectivity be assessed per-epoch?

**Recommendation:** Extract features separately per epoch, let clustering discover which epochs matter most.

### Q3: Handling the "Almost Classic" Neurons

**Issue:** User reports some neurons are "clearly classic" but get classified as interneuron due to weak ipsilateral responses.

**Question:** How to handle neurons where:
- Contralateral response is 10x ipsilateral response
- But ipsilateral response is still significant

**Potential approaches:**
1. **Magnitude threshold:** Require ipsi response to be >X% of contra to count as "opposing"
2. **Selectivity index:** Use continuous measure rather than binary
3. **Clustering:** Let the data reveal whether these form a distinct group

### Q4: Handling Suppression

**Issue:** Current approach flags ANY contralateral suppression as "interneuron-like."

**Question:** Should suppression be weighted by:
- Magnitude relative to activation?
- Epoch (visual vs saccade)?
- Location (best vs non-best)?

**Recommendation:** Include suppression strength as a continuous feature, let clustering reveal whether suppression-dominant neurons form a distinct group.

### Q5: Relationship to 4-Factors Selectivity

**Issue:** The ultimate goal is to understand how different SC neuron types encode the task factors.

**Question:** Should factor selectivity influence classification, or should classification be factor-blind?

**Recommendation:** Classification should be factor-blind (based only on spatial/temporal response patterns). Factor selectivity is then analyzed as a dependent variable to test hypotheses about functional roles.

---

## 6. Alternative Approaches to Consider

### Alternative 1: Rule-Based Refinement

Instead of clustering, refine the current rule-based approach:

```matlab
% Spatial selectivity index
SSI = (contra_activation - ipsi_activation) / (contra_activation + ipsi_activation);

% Classifications
if SSI > 0.5 && no_suppression
    class = 'classic_selective';
elseif SSI > 0 && SSI <= 0.5
    class = 'classic_broad';
elseif has_suppression && SSI > 0
    class = 'suppressed_selective';
elseif SSI <= 0
    class = 'ipsi_preferring';
end
```

**Pros:** Transparent, easy to interpret
**Cons:** May not capture complex patterns, requires tuning thresholds

### Alternative 2: Supervised Classification

If expert annotations are available:

1. User annotates subset of neurons
2. Train classifier (SVM, random forest) on annotated examples
3. Apply to remaining neurons
4. Iterate with active learning

**Pros:** Directly captures expert intuition
**Cons:** Requires significant annotation effort, may not generalize

### Alternative 3: Hybrid Approach

1. Use clustering to identify candidate groupings
2. Present cluster exemplars to user for interpretation
3. Refine cluster boundaries based on user feedback
4. Codify final rules for reproducibility

**This is the recommended approach.**

---

## 7. Success Criteria

The new classification will be considered successful if:

1. **Reproducibility:** Same neurons classified consistently across runs (100%)

2. **Face validity:** Expert (user) agrees with classification for >80% of neurons

3. **Separability:** Clusters are statistically distinct (silhouette score > 0.4)

4. **Biological interpretability:** Each cluster has a plausible biological interpretation

5. **Functional relevance:** Classification improves understanding of factor encoding
   - Different classes show different patterns of 4-factors selectivity
   - Classification explains variance not captured by current approach

6. **Parsimony:** Final number of categories is manageable (3-7, not 15+)

---

## 8. Timeline and Milestones

| Week | Milestone | Deliverables |
|------|-----------|--------------|
| 1 | Feature extraction | `extract_sc_response_features.m`, exploratory figures |
| 1-2 | Clustering analysis | Optimal k determination, cluster assignments |
| 2 | Interpretation | Cluster labels, representative neuron examples |
| 2-3 | Validation | Expert agreement rate, cross-validation results |
| 3 | Integration | Updated `screen_sc_neurons.m`, documentation |

---

## 9. Files to Create/Modify

### New Files

| File | Purpose |
|------|---------|
| `code/extract_sc_response_features.m` | Feature extraction from cached results |
| `code/cluster_sc_neurons.m` | Unsupervised clustering pipeline |
| `code/visualize_sc_clusters.m` | Cluster visualization and interpretation |
| `code/validate_sc_classification.m` | Validation against expert annotation |
| `docs/sc_neuron_cluster_interpretation.md` | Documentation of final categories |

### Modified Files

| File | Changes |
|------|---------|
| `code/screen_sc_neurons.m` | Add new classification, maintain backward compatibility |
| `code/aggregate_analysis_results.m` | Aggregate new classification fields |
| `docs/design/analysis_rationale.md` | Update neuron classification section |

---

## 10. Open Questions for Discussion

Before implementation, the following should be discussed:

1. **Priority of interpretability vs accuracy:** If clustering suggests 8 distinct groups, do we keep all 8 or collapse to 4-5 for simplicity?

2. **Treatment of excluded neurons:** Should we attempt to classify neurons that fail basic criteria (low FR, sparse), or focus only on "includable" neurons?

3. **Integration with existing analyses:** Should we re-run all downstream analyses with the new classification, or treat it as supplementary?

4. **Manual override mechanism:** Should we maintain the ability to manually override algorithmic classifications for specific neurons?

5. **Version control:** How do we handle the transition? Keep both classifications in parallel?

---

## 11. References

1. Phongphanphanee, P., et al. (2023). Functional classes in the superior colliculus. [PubMed 37202508](https://pubmed.ncbi.nlm.nih.gov/37202508/)

2. Lee, K.H., et al. (2021). Inhibitory neurons in the superior colliculus mediate selection of spatially-directed movements. [Nature Communications Biology](https://www.nature.com/articles/s42003-021-02248-1)

3. Whyland, K.L., et al. (2023). Parvalbumin and GABA Microcircuits in the Mouse Superior Colliculus. [Frontiers in Neural Circuits](https://www.frontiersin.org/articles/10.3389/fncir.2018.00035)

4. Drion, G., et al. (2015). Towards automatic classification of neurons. [PLOS Computational Biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC4417416/)

5. Goetz, L., et al. (2022). Functional cell types in the mouse superior colliculus. [eLife](https://elifesciences.org/articles/82367)

6. Baden, T., et al. (2016). The functional diversity of retinal ganglion cells in the mouse. [Nature](https://www.nature.com/articles/nature16468) - Example of clustering-based cell type identification

---

## Appendix A: Current Feature Availability

From `cached_friedman_results` in `screen_sc_neurons.m`:

```matlab
cached_friedman_results = struct(...
    'cluster_ids',           [],        % [nClusters × 1]
    'sig_by_location',       [],        % [nClusters × nGroups × 3] (vis, delay, sac)
    'suppression_by_location', [],      % [nClusters × nGroups × 3]
    'group_mean_frs',        {},        % {1 × nGroups} cell, each [nClusters × 4] (base, vis, delay, sac)
    'location_angles',       [],        % [nGroups × 1] angles in degrees
    'epoch_frs',             [],        % [nClusters × 4 × nTrials] raw firing rates
    'computation_date',      ''
);
```

This provides:
- Significance of activation/suppression at each location × epoch
- Mean firing rates by location × epoch
- Trial-by-trial firing rates for variance estimation

**Missing but potentially useful:**
- Spike waveform features (duration, asymmetry)
- ISI statistics (CV, burst index)
- Response latency estimates
