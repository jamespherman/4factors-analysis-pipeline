# Analysis Rationale Document: 4factors Analysis Pipeline

This document connects analyses to scientific hypotheses, documents design rationale, and serves as a planning blueprint for the 4factors analysis pipeline. It is intended to explain the "why" behind each analysis and how outcomes would be interpreted.

---

## Table of Contents

1. [Scientific Context and Hypotheses](#1-scientific-context-and-hypotheses)
2. [Analysis Architecture Overview](#2-analysis-architecture-overview)
3. [Neuron Classification and Screening](#3-neuron-classification-and-screening)
4. [Single-Neuron Analyses](#4-single-neuron-analyses)
5. [Population-Level Analyses](#5-population-level-analyses)
6. [SNc Subregion Analysis](#6-snc-subregion-analysis)
7. [SC-SNc Interaction Analysis](#7-sc-snc-interaction-analysis)
8. [Behavioral Analyses](#8-behavioral-analyses)
9. [Visualization and Verification](#9-visualization-and-verification)
10. [Implementation Priorities](#10-implementation-priorities)
11. [Appendix: Technical Specifications](#appendix-technical-specifications)

---

## 1. Scientific Context and Hypotheses

### 1.1 Background

The Superior Colliculus (SC) integrates multiple sources of information to generate a "priority map" that guides saccadic eye movements. These sources include:

- **Stimulus-driven factors**: Physical salience, object identity (face vs non-face vs geometric)
- **Goal-directed factors**: Expected reward value, spatial probability (learned expectation)

The Substantia Nigra pars compacta (SNc) provides dopaminergic input to the striatum and is implicated in reward prediction and learning. The SNc is not homogeneous—different subregions project to different parts of the caudate nucleus:

- **Rostral ventromedial SNc (rvmSNc)**: Projects to the head of the caudate, which receives cortical inputs related to goal-directed behavior
- **Caudal dorsolateral SNc (cdlSNc)**: Projects to the tail of the caudate, which receives cortical inputs related to visual processing and stimulus-driven responses

### 1.2 Primary Hypotheses

**Hypothesis 1: Factor Separability in SC**
> Stimulus-driven factors (salience, identity) and goal-directed factors (reward, probability) are encoded in separable subspaces within SC population activity. That is, a decoder trained to discriminate high vs low reward cannot generalize to discriminate high vs low salience, and vice versa.

**Hypothesis 2: SNc Subregion Specialization**
> SNc dopamine neurons show subregion-specific factor encoding:
> - rvmSNc neurons preferentially encode goal-directed factors (reward, probability)
> - cdlSNc neurons preferentially encode stimulus-driven factors (salience, identity)

**Hypothesis 3: SC→SNc Information Routing**
> SC transmits factor-specific information to SNc in a pathway-selective manner:
> - Stimulus-driven information flows preferentially to cdlSNc
> - Goal-directed information flows preferentially to rvmSNc

**Hypothesis 4: SC Interneuron Non-Selectivity**
> Putative SC interneurons (identified by visual suppression rather than excitation) are non-selective for task factors, consistent with providing broad inhibition rather than feature-specific processing.

### 1.3 The 4factors Task Design

The gSac_4factors task is a memory-guided saccade paradigm with a 2×2×2×2 factorial design (though not fully crossed due to the salience/identity confound):

| Factor | Levels | Manipulation |
|--------|--------|--------------|
| **Stimulus Type** | Face, Non-face, Bullseye | Different visual targets |
| **Salience** | High, Low | Background contrast (bullseye only) |
| **Reward** | High, Low | Juice duration, spatially mapped |
| **Probability** | High, Low | Block-contingent spatial expectation |

**Critical design feature**: Salience only varies for bullseye trials (images are presented on isoluminant gray). This necessitates a **two-model approach**:
- **Image Trials Model**: Reward × Probability × Identity
- **Bullseye Trials Model**: Reward × Probability × Salience

---

## 2. Analysis Architecture Overview

### 2.1 Pipeline Stages

```
┌─────────────────────────────────────────────────────────────────┐
│                    STAGE 1: Per-Session                         │
├─────────────────────────────────────────────────────────────────┤
│  1. Load session_data.mat (from preprocessing pipeline)         │
│  2. Neuron screening and classification                         │
│  3. Core data preparation (spike binning, trial filtering)      │
│  4. Single-neuron analyses (ROC, baseline comparison)           │
│  5. Population analyses (ANOVA, decoding)                       │
│  6. Generate per-neuron diagnostic PDFs                         │
│  7. Save results to session_data.analysis                       │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    STAGE 2: Aggregation                         │
├─────────────────────────────────────────────────────────────────┤
│  1. Load all processed session_data files                       │
│  2. Concatenate results across sessions                         │
│  3. Compute population-level statistics                         │
│  4. Save aggregated_analysis_data.mat                           │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    STAGE 3: Plotting                            │
├─────────────────────────────────────────────────────────────────┤
│  1. Load aggregated data                                        │
│  2. Generate population summary figures                         │
│  3. Generate scatter plots with neuron identification           │
│  4. Save figures to /figures/ directory                         │
└─────────────────────────────────────────────────────────────────┘
```

### 2.2 Data Flow

```
session_data.mat (from preprocessing)
    │
    ├── trialInfo (table): Trial parameters
    ├── eventTimes (struct): Timestamps for key events
    ├── spikes (struct): times, clusters, cluster_info, waveforms
    └── metadata (struct): session identifiers, grid coordinates
                              │
                              ▼
session_data.analysis (added by this pipeline)
    │
    ├── selected_neurons: Logical vector of included neurons
    ├── neuron_classification: Categorical type for each neuron
    ├── core_data: Binned spike rates aligned to events
    ├── roc_comparison: Time-resolved ROC results
    ├── window_roc: Epoch-based ROC values
    ├── baseline_comparison: Post-event vs baseline ROC
    ├── anova_results: Time-resolved ANOVA p-values
    ├── behavioral_results: Behavioral ANOVA results
    └── population_decoding: SVM accuracy and generalization
```

### 2.3 Condition Definitions

All condition masks are defined in `define_task_conditions.m`, which serves as the single source of truth. Key masks include:

**Spatial:**
- `is_contralateral_target`: Target in the hemifield contralateral to recording site
- `is_ipsilateral_target`: Target in ipsilateral hemifield
- `is_opposite_rf`: Target at the location opposite to the primary RF location

**Factorial:**
- `is_high_reward` / `is_low_reward`: Based on reward duration threshold (200ms)
- `is_high_salience` / `is_low_salience`: Based on salience parameter (bullseye only)
- `is_face_target` / `is_nonface_target`: Based on stimulus type
- `is_high_probability` / `is_low_probability`: Based on block-contingent design

**Trial Type:**
- `is_image_target`: Face or non-face trials
- `is_bullseye_target`: Bullseye trials

**SNc Subregion:**
- `is_rvmSNc`: Grid Y coordinate is 3 or 4
- `is_cdlSNc`: Grid Y coordinate is 1 or 2

---

## 3. Neuron Classification and Screening

### 3.1 Current Approach

**SC Neurons** (`screen_sc_neurons.m`):
- Uses Friedman test across four epochs (baseline, visual, delay, saccade)
- If significant, performs Bonferroni-corrected Wilcoxon tests comparing baseline to each other epoch
- Selects neurons showing significant **increase** in any epoch with max FR > 5 sp/s

**DA Neurons** (`screen_da_neurons.m`):
- Selects based on: baseline FR < 20 sp/s AND waveform duration > 0.6 ms

### 3.2 Proposed Expansion: Neuron Classification

**Rationale**: The current SC screening excludes neurons that are suppressed by visual stimuli or unmodulated. Hypothesis 4 specifically concerns these excluded neurons.

**Proposed Classification Scheme for SC:**

| Category | Criteria | Interpretation |
|----------|----------|----------------|
| `visually_excited` | Significant increase in visual epoch vs baseline (contralateral targets) | Classical SC visual neurons |
| `visually_suppressed` | Significant decrease in visual epoch vs baseline (contralateral targets) | Putative interneurons |
| `delay_modulated` | Significant change in delay epoch (either direction) | Memory-related activity |
| `saccade_modulated` | Significant change in saccade epoch (either direction) | Motor-related activity |
| `unmodulated` | No significant modulation in any epoch | Non-task-related or tonically active |

**Critical Note on Suppression**: The key criterion for `visually_suppressed` is suppression for **contralateral** target presentations. Classical visually-excited SC neurons are excited by contralateral targets and may be suppressed by ipsilateral targets—this is expected spatial selectivity, not interneuron-like behavior. True putative interneurons show suppression when targets appear in what would be the excitatory RF for projection neurons (i.e., contralateral hemifield).

**Statistical Threshold**: Same as for excitation—Bonferroni-corrected Wilcoxon rank-sum test (alpha = 0.05/3 for 3 comparisons). No minimum magnitude requirement; the key is statistical reliability of the suppression.

**Implementation**: Replace the single `selected_neurons` logical vector with:
- `neuron_type`: Categorical variable with the above categories
- `include_in_factor_analysis`: Logical vector (true for visually_excited, visually_suppressed, delay_modulated, saccade_modulated)

**Hypothesis 4 Test**: For neurons classified as `visually_suppressed`:
1. Compute factor-specific ROC values (same as for excited neurons)
2. Test whether ROC values differ from 0.5 (chance)
3. Prediction: Suppressed neurons will show ROC ≈ 0.5 for all factors (non-selective)

### 3.3 SNc Subregion Labeling

**Rule** (from manifest):
- First, verify `brain_area == 'SNc'` (SC neurons are not labeled with subregion)
- Parse Y coordinate from `grid_hole` column (format: `(X, Y)`)
- Y ∈ {3, 4} → `rvmSNc`
- Y ∈ {1, 2} → `cdlSNc`

**Implementation**: Add `snc_subregion` field to `session_data.metadata` during session loading (set to empty string for SC sessions), and create corresponding condition masks in `define_task_conditions.m`.

---

## 4. Single-Neuron Analyses

### 4.1 Time-Resolved ROC (Existing)

**Purpose**: Track moment-by-moment discriminability between conditions.

**Method**: For each neuron and each time bin, compute ROC AUC comparing spike counts between two conditions (e.g., high vs low reward). Significance determined by permutation test (200 shuffles).

**Output**: `[n_neurons × n_bins]` matrices for AUC and significance.

**Interpretation**: AUC > 0.5 indicates higher firing for condition 1; AUC < 0.5 indicates higher firing for condition 2. Significance indicates reliable discrimination.

### 4.2 Window-Based ROC

**Purpose**: Provide a single summary metric per neuron per factor per epoch, suitable for scatter plots and cross-neuron comparisons.

**Rationale for window-based over max-in-window**: Using max AUC within a window can be misleading when comparing across factors. If neuron A has max AUC for reward at 100ms post-target but max AUC for salience at 300ms post-target, comparing these values conflates timing with magnitude. A fixed window ensures we're comparing the same temporal epoch.

**Proposed Windows**:

| Epoch | Event | SC Window (ms) | SNc Window (ms) |
|-------|-------|----------------|-----------------|
| Visual | targetOn | [50, 250] | [100, 300] |
| Delay | fixOff | [-200, 0] | [-200, 0] |
| Peri-saccadic | saccadeOnset | [-50, 50] | [-50, 50] |
| Post-reward | reward | [100, 400] | [100, 400] |

**Rationale for area-specific windows**: SNc visual responses are slower than SC (~100ms vs ~50ms latency), so the visual window is shifted later for SNc neurons.

**Method**: 
1. For each neuron, compute mean spike count in the window for each trial
2. Compute ROC AUC comparing the two conditions
3. Compute significance via permutation test
4. Store: AUC value, p-value, and 95% CI

**Output Structure**:
```matlab
window_roc.visual.reward.auc      % [n_neurons × 1]
window_roc.visual.reward.p        % [n_neurons × 1]
window_roc.visual.reward.ci       % [n_neurons × 2]
window_roc.visual.salience.auc
... etc for all epoch × factor combinations
```

**Interpretation for Hypothesis 2**:

The key scatter plots compare stimulus-driven factors against goal-directed factors:

**Primary Comparisons (Cross-Category)**:
| Plot | X-axis | Y-axis | Hypothesis 2 Prediction |
|------|--------|--------|-------------------------|
| 1 | Reward ROC | Salience ROC | rvmSNc clusters high-X; cdlSNc clusters high-Y |
| 2 | Reward ROC | Identity ROC | rvmSNc clusters high-X; cdlSNc clusters high-Y |
| 3 | Probability ROC | Salience ROC | rvmSNc clusters high-X; cdlSNc clusters high-Y |
| 4 | Probability ROC | Identity ROC | rvmSNc clusters high-X; cdlSNc clusters high-Y |

**Secondary Comparisons (Within-Category)**:
| Plot | X-axis | Y-axis | Expected Pattern |
|------|--------|--------|------------------|
| 5 | Reward ROC | Probability ROC | Positive correlation (both goal-directed) |
| 6 | Salience ROC | Identity ROC | Positive correlation (both stimulus-driven) |

Each scatter plot should be generated for at least the Visual epoch. Saccade epoch plots may also be informative.

**Summary Prediction**: If Hypothesis 2 is correct, the primary cross-category plots should show **different clustering patterns** for rvmSNc vs cdlSNc neurons, while the within-category plots should show **similar positive correlations** for both subregions.

### 4.3 Baseline Comparison (Existing)

**Purpose**: Test whether neurons are modulated by the task at all (activity vs baseline).

**Method**: ROC comparing post-event firing to pre-event baseline for each condition separately.

**Interpretation**: Identifies neurons that increase (or decrease) firing in response to task events, independent of factor selectivity.

---

## 5. Population-Level Analyses

### 5.1 N-way ANOVA (Existing)

**Purpose**: Test main effects and interactions of task factors on neural activity.

**Two-Model Approach**:
- **Image Trials**: `anovan(firing_rate, {reward, probability, identity})`
- **Bullseye Trials**: `anovan(firing_rate, {reward, probability, salience})`

**Output**: Time-resolved p-values and F-statistics for each term (main effects + interactions).

**Interpretation for Hypothesis 1**:
- Significant main effects indicate factor encoding
- Significant interactions would suggest factors are not fully independent
- Prediction: Main effects present; interactions weak or absent

### 5.2 Population Decoding (Existing)

**Purpose**: Test whether factor information is present at the population level and whether it generalizes across factors/epochs.

**Method**: Linear SVM with 10-fold cross-validation.

**Test Types**:

| Test Type | Training | Testing | Question |
|-----------|----------|---------|----------|
| `standard` | Factor X, Epoch A | Same (cross-validated) | Can we decode factor X? |
| `cross_time` | Factor X, Epoch A | Factor X, Epoch B | Does encoding persist across epochs? |
| `cross_factor` | Factor X | Factor Y | Does factor X encoding generalize to factor Y? |

**Interpretation for Hypothesis 1**:
- High `standard` accuracy: Factor is encoded
- Low `cross_factor` accuracy: Factors are encoded independently (separable)
- Prediction: `standard` accuracy >> `cross_factor` accuracy for stimulus-driven vs goal-directed factors

---

## 6. SNc Subregion Analysis

### 6.1 Labeling

**Source**: `grid_hole` column in session manifest, format `(X, Y)`

**Rule**:
```matlab
coords = sscanf(grid_hole_str, '(%f, %f)');
grid_y = coords(2);
if grid_y >= 3
    snc_subregion = 'rvmSNc';
else
    snc_subregion = 'cdlSNc';
end
```

### 6.2 Hypothesis 2 Tests

**Primary Analysis**: Compare factor selectivity between subregions.

**Metrics to compare**:
- Mean window-based ROC AUC for each factor
- Proportion of neurons with significant selectivity for each factor
- Population decoding accuracy for each factor

**Statistical Tests**:
- Wilcoxon rank-sum test comparing ROC distributions between subregions
- Chi-square test for independence of factor selectivity: Tests whether the proportion of neurons selective for **both** factors in a pair (e.g., salience AND reward) differs from what would be expected under independence (i.e., P(selective for both) vs P(selective for A) × P(selective for B)). A significant result indicates that selectivity for the two factors is not independent across neurons.
- Bootstrap comparison of decoding accuracies

**Visualization**:
- Scatter plot: Reward ROC vs Salience ROC, colored by subregion
- Bar plot: Mean ROC by factor × subregion
- Violin plot: ROC distributions by subregion

**Predictions**:
- rvmSNc: Higher ROC for reward and probability; lower ROC for salience and identity
- cdlSNc: Higher ROC for salience and identity; lower ROC for reward and probability

---

## 7. SC-SNc Interaction Analysis

### 7.1 Overview

**Goal**: Test whether SC transmits factor-specific information to SNc subregions selectively.

**Data Requirement**: Simultaneous SC + SNc recordings (available: 4 sessions in Feynman)

### 7.2 Session Organization

Simultaneous SC and SNc recordings are stored as **separate `_session_data.mat` files** linked by the `session_group_id` column in the manifest. Sessions sharing the same `session_group_id` and `raw_filename_base` were recorded simultaneously.

**Available Simultaneous Sessions:**

| Date | Group ID | SNc Session | SC Session | SNc Grid Y | SNc Subregion |
|------|----------|-------------|------------|------------|---------------|
| 08/12/2025 | Feynman_08_12_2025 | Feynman_08_12_2025_SNc | Feynman_08_12_2025_SC | 3 | rvmSNc |
| 08/13/2025 | Feynman_08_13_2025 | Feynman_08_13_2025_SNc | Feynman_08_13_2025_SC | 1 | cdlSNc |
| 08/14/2025 | Feynman_08_14_2025 | Feynman_08_14_2025_SNc | Feynman_08_14_2025_SC | 3 | rvmSNc |
| 08/15/2025 | Feynman_08_15_2025 | Feynman_08_15_2025_SNc | Feynman_08_15_2025_SC | 2 | cdlSNc |

**Balance**: 2 sessions with rvmSNc, 2 sessions with cdlSNc—well-balanced for subregion comparisons.

### 7.3 Approach 1: Correlation of Factor Selectivity (Simpler)

**Method**:
1. For each simultaneous session, compute:
   - SC population decoding accuracy for each factor
   - SNc subregion-specific mean ROC for each factor
2. Across the 4 sessions, correlate SC factor decoding with SNc factor encoding
3. Test whether correlations differ by SNc subregion

**Predictions**:
- SC salience encoding correlates with cdlSNc salience encoding (but not rvmSNc)
- SC reward encoding correlates with rvmSNc reward encoding (but not cdlSNc)

**Limitation**: Only 4 data points per correlation. Results will be suggestive, not definitive.

### 7.4 Approach 2: Canonical Correlation Analysis (More Sophisticated)

**Method**:
1. For each trial, construct:
   - SC population vector: Firing rates of all SC neurons in a given epoch
   - SNc population vector: Firing rates of all SNc neurons in the same epoch
2. Perform CCA to find linear combinations that maximize correlation
3. Examine the canonical correlations and weights

**Analysis Extensions**:
a. **CCA weight vs factor selectivity**: Do SNc neurons with high CCA weights (strongly coupled to SC) tend to have specific factor selectivity?
b. **Subregion-specific CCA**: Perform CCA separately for rvmSNc and cdlSNc. Do canonical correlations differ?
c. **Factor-conditioned CCA**: Perform CCA on residuals after regressing out factor effects. Does SC-SNc coupling exist beyond shared factor encoding?

**Visualization (Your Proposed Analysis)**:
- Scatter plot: Mean SNc goal-directed ROC vs CCA weight, colored by subregion
- Scatter plot: Mean SNc stimulus-driven ROC vs CCA weight, colored by subregion

**Predictions** (testing pathway-selective information routing):

| Subregion | Factor Type | ROC vs CCA Weight Correlation | Interpretation |
|-----------|-------------|------------------------------|----------------|
| rvmSNc | Goal-directed | **Significant positive** | SC transmits goal-directed info to rvmSNc |
| rvmSNc | Stimulus-driven | **Not significant** | SC does NOT transmit stimulus-driven info to rvmSNc |
| cdlSNc | Stimulus-driven | **Significant positive** | SC transmits stimulus-driven info to cdlSNc |
| cdlSNc | Goal-directed | **Not significant** | SC does NOT transmit goal-directed info to cdlSNc |

This double-dissociation pattern would support Hypothesis 3 (pathway-selective information routing).

**Interpretation Challenges**:
- CCA finds shared variance, not causal direction
- Shared inputs (e.g., from cortex) could drive SC-SNc correlations without direct SC→SNc transmission
- Consider using time-lagged CCA or Granger-style analyses if temporal resolution permits

### 7.5 Implementation Considerations

**Linking Simultaneous Sessions**:
```matlab
% Find paired sessions by matching session_group_id
manifest = readtable('session_manifest.csv');
unique_groups = unique(manifest.session_group_id);

for i = 1:length(unique_groups)
    group_id = unique_groups{i};
    group_rows = strcmp(manifest.session_group_id, group_id);
    
    if sum(group_rows) == 2  % Paired SC + SNc
        sessions = manifest(group_rows, :);
        sc_session = sessions(strcmp(sessions.brain_area, 'SC'), :);
        snc_session = sessions(strcmp(sessions.brain_area, 'SNc'), :);
        % Process paired sessions...
    end
end
```

---

## 8. Behavioral Analyses

### 8.1 Current Analyses

**Reaction Time (Saccadic Reaction Time, SRT)**:
- DV: `eventTimes.saccadeOnset - eventTimes.fixOff`
- Model: N-way ANOVA with factors from the two-model approach

**Peak Velocity**:
- DV: `trialInfo.peakVel`
- Model: Same as RT

**Endpoint Error**:
- DV: Distance from saccade endpoint to target
- Model: Same as RT

### 8.2 Additional Analyses

**Error Rate**:
- DV: Proportion of trials not resulting in reward, by condition
- Method: Logistic regression or chi-square test
- Interpretation: Higher error rates may indicate increased difficulty or reduced motivation

**Main Sequence Analysis**:
- Standard approach: `PeakVel = a × Amplitude^b`
- In this task, all targets are presented at equal eccentricity, so amplitude is approximately constant across conditions (varying only due to landing accuracy)
- This simplifies the analysis: peak velocity can be analyzed directly without amplitude confound

**Proposed Method**:
1. Analyze peak velocity directly with condition factors (reward, probability, etc.)
2. Optionally include saccade amplitude as a covariate to account for trial-to-trial landing variability
3. Interpretation: Higher peak velocity = more vigorous saccade = potentially higher motivation/value

### 8.3 Methodological Consideration: Saccade Direction Effects

**Concern**: Oculomotor metrics (SRT, peak velocity, accuracy) may differ systematically between leftward and rightward saccades due to biomechanical or neural asymmetries. If large, these direction effects could mask factor-related effects.

**Recommended Approaches**:

1. **Include direction as a factor**: Add saccade direction (left/right) as an additional factor in the behavioral ANOVA. This allows estimation of direction main effects and direction × factor interactions.

2. **Within-direction analysis**: If direction effects are very large, analyze leftward and rightward saccades separately, then compare effect sizes for task factors within each direction.

3. **Direction as covariate**: Include direction as a covariate to partial out its effects when testing factor effects.

**Implementation**: The first approach (direction as factor) is preferred as it provides the most complete picture and allows detection of direction × factor interactions that could be scientifically meaningful.

---

## 9. Visualization and Verification

### 9.1 Design Philosophy

Every population-level metric should have a corresponding visualization in the single-neuron summaries. This enables:
- **Verification**: Confirm that population effects are consistent with single-neuron data
- **Example Selection**: Identify neurons that exemplify population trends for figures
- **Anomaly Detection**: Spot neurons that deviate from expected patterns

### 9.2 Existing Population-Level Figures

#### 9.2.1 Aggregated ANOVA Plot (`plot_aggregated_anova.m`)

**Purpose**: Visualize proportion of neurons with significant main effects and interactions over time.

**Layout**:
- Rows: ANOVA terms (main effects and interactions from both models)
  - Image Trials Model: reward, probability, identity, and their interactions
  - Bullseye Trials Model: reward, probability, salience, and their interactions
- Columns: Alignment events (fixOn, targetOn, fixOff, saccadeOnset, reward)

**Plot Type**: Line plot showing mean proportion significant across sessions (bold) with individual session traces (faint).

**Current Status**: Implemented, functional.

**Potential Improvements**: 
- Could add shaded regions for SEM instead of individual traces
- Could separate Image and Bullseye models into different figures or add model indicator

---

#### 9.2.2 Aggregated Baseline Comparison Plot (`plot_aggregated_baseline_comparison.m`)

**Purpose**: Visualize proportion of neurons with significant increase/decrease from baseline, by condition.

**Layout**:
- Rows: Conditions (is_high_reward, is_low_reward, is_high_salience, etc.)
- Columns: Alignment events

**Plot Type**: Bidirectional `barStairsFill` (positive = increase, negative = decrease)

**Current Status**: Implemented, functional.

**Potential Improvements**:
- Group rows by factor type (all reward conditions together, etc.)

---

#### 9.2.3 Aggregated ROC Comparison Plot (`plot_aggregated_roc_comparison.m`)

**Purpose**: Visualize proportion of neurons preferring one condition over another, over time.

**Layout**:
- Rows: Factor comparisons (reward, salience, identity, probability)
- Columns: Alignment events

**Plot Type**: Bidirectional `barStairsFill` (positive = prefers cond2, negative = prefers cond1)

**Current Status**: Implemented, functional.

**Interpretation for Hypothesis Testing**:
- Large proportions indicate population-level factor encoding
- Timing of selectivity reveals when information emerges (visual vs delay vs peri-saccadic)

---

#### 9.2.4 Aggregated Behavioral Plot (`plot_aggregated_behavior.m`)

**Purpose**: Visualize proportion of sessions with significant behavioral effects.

**Layout**:
- Rows: Behavioral measures (RT, peak velocity, endpoint error)
- Columns: Main effects vs interactions

**Plot Type**: Grouped bar plot (Image model vs Bullseye model)

**Current Status**: Implemented, functional.

**Potential Improvements**:
- Add actual effect sizes (not just proportion significant)
- Add within-session distributions

---

#### 9.2.5 Aggregated Decoding Plot (`plot_aggregated_decoding.m`)

**Purpose**: Visualize decoder performance and generalization.

**Current Layout**: Flat grid (4 columns × N rows) with one subplot per generalization test.

**Plot Type**: Scatter plot (X = standard CV accuracy, Y = generalization accuracy), one point per session.

**Current Status**: Implemented but **suboptimally organized**.

**Issues Identified**:

1. **Lack of Logical Grouping**: Cross-factor and cross-time tests are intermixed in the grid. Tests should be grouped by type to facilitate interpretation.

2. **Missing Organization by Factor**: Within cross-factor tests, there's no clear visual organization showing "train on X, test on Y" relationships. A matrix layout would be more intuitive.

3. **No Summary Statistics**: The plot shows raw data but doesn't highlight the key finding (e.g., whether cross-factor generalization is systematically worse than within-factor).

4. **Dense Layout**: With 4×2 cross-factor combinations × 2 epochs + 4×2 cross-time combinations = 24+ tests, the grid becomes crowded.

**Reorganization**:

The decoding visualization is split into two figures for clarity:

**Figure 1: Cross-Factor Generalization** (matrix layout)
- Rows: Training factor (Reward, Salience, Identity, Probability)
- Columns: Testing factor
- Each cell shows Visual and Saccade epoch results
- Diagonal cells = standard CV (within-factor baseline)
- Off-diagonal cells = cross-factor generalization

**Figure 2: Cross-Time Generalization**
- Rows: Factors (Reward, Salience, Identity, Probability)
- Columns: Direction (Visual→Saccade, Saccade→Visual)
- Includes standard CV for reference

This organization makes the key question immediately visible: Do off-diagonal cells (cross-factor) show systematically lower accuracy than diagonal cells (within-factor)?

---

### 9.3 Population-Level Figures: Factor Selectivity

#### 9.3.1 ROC Scatter Plots

**Purpose**: Compare single-neuron selectivity across factors, enabling identification of neurons for follow-up.

**Layout**: Multi-panel grid showing all pairwise factor comparisons:
- 6 panels for 4 factors: Reward×Salience, Reward×Probability, Reward×Identity, Salience×Probability, Salience×Identity, Probability×Identity
- Each panel: X = Factor A ROC, Y = Factor B ROC
- Points colored by brain area (SC vs SNc) or SNc subregion (rvmSNc vs cdlSNc)

**Plot Type**: Scatter plot with hover tooltips for neuron identification.

**Data Requirements**:
- Window-based ROC values per neuron (from `analyze_window_roc.m`)
- Neuron metadata (session_id, cluster_id, brain_area, snc_subregion)

**Tooltip Content**:
```
Session: Feynman_08_12_2025_SNc
Cluster: 14
Type: visually_excited
Subregion: rvmSNc
Reward ROC: 0.72
Salience ROC: 0.48
```

**Interpretation for Hypothesis 2**:
- rvmSNc neurons should cluster toward high Reward/Probability ROC
- cdlSNc neurons should cluster toward high Salience/Identity ROC
- Separation between subregions supports the hypothesis

**CSV Backup**:
```
/figures/scatter_data/roc_reward_vs_salience_visual.csv
Columns: session_id, cluster_id, brain_area, snc_subregion, neuron_type, reward_roc, salience_roc
```

---

#### 9.3.2 SNc Subregion Comparison Figures

**Purpose**: Directly compare factor selectivity between rvmSNc and cdlSNc.

**Layout**:
- **Panel A**: Grouped bar plot showing mean ROC by factor × subregion
- **Panel B**: Violin plots showing full ROC distribution by subregion
- **Panel C**: Scatter plot with subregion color coding

**Statistical Annotations**: p-values from Wilcoxon rank-sum tests comparing subregions.

---

#### 9.3.3 SC-SNc Interaction Figures (Simultaneous Sessions)

**Purpose**: Visualize relationship between SC population encoding and SNc factor selectivity.

**Figure 1: Factor Selectivity Correlation**
- Scatter plot: SC population decoding accuracy vs SNc mean ROC
- Separate panels for stimulus-driven vs goal-directed factors
- Color by SNc subregion
- Note: Only 4 data points (sessions)

**Figure 2: CCA Weight Analysis (if implemented)**
- Scatter plot: CCA weight vs factor ROC for each SNc neuron
- Color by subregion
- Panel A: Stimulus-driven factors
- Panel B: Goal-directed factors

---

### 9.4 Single-Neuron Summary PDFs

**Guiding Principle**: For any population scatter plot, the user should be able to:
1. Identify a neuron of interest via hover tooltip
2. Open that neuron's PDF
3. Visually confirm the underlying data matches the summary statistics

**Scope**: Generate PDFs for all modulated and suppressed neurons (not unmodulated).

**Single-Page Layout** (due to MATLAB multi-page PDF limitations):

The layout uses a 5-column structure to present all relevant information on one page:

| Column 1 | Column 2 | Column 3 | Column 4 | Column 5 |
|----------|----------|----------|----------|----------|
| Spatial: targetOn (jph) | Spatial: saccadeOnset (jph) | Factors: targetOn | Factors: saccadeOnset | Factors: reward |

**Row Organization** (top to bottom):

**Header Row**: Neuron metadata spanning full width
- Cluster ID, neuron type, brain area, SNc subregion (if applicable)
- Baseline FR, waveform duration, Phy quality label

**Rows 1-2**: Overview panels (Columns 1-2) | Factor PSTHs: Reward (Columns 3-5)
- Col 1-2: Waveform (multi-channel heatmap or mean waveform)
- Col 1-2: ISI histogram with refractory violation %
- Col 3-5: Reward high vs low PSTHs for each alignment event

**Rows 3-4**: Spatial tuning (Columns 1-2) | Factor PSTHs: Salience (Columns 3-5)
- Col 1-2: Raster + PSTH for contra vs ipsi (gSac_jph trials)
- Col 3-5: Salience high vs low PSTHs (bullseye trials only)

**Rows 5-6**: Spatial tuning continued (Columns 1-2) | Factor PSTHs: Probability (Columns 3-5)
- Col 1-2: Raster + PSTH for second spatial alignment
- Col 3-5: Probability high vs low PSTHs

**Rows 7-8**: Summary statistics (Columns 1-2) | Factor PSTHs: Identity (Columns 3-5)
- Col 1-2: ROC summary table (factors × epochs)
- Col 3-5: Identity face vs non-face PSTHs (image trials only)

**Color Scheme**:
- Spatial tuning: Blue = Contralateral, Red = Ipsilateral
- Factor comparisons: Blue = High condition, Red = Low condition
- Shaded regions = SEM

**ROC Summary Table** (embedded in layout):

| Factor | Visual | Delay | Saccade | Reward |
|--------|--------|-------|---------|--------|
| Reward | 0.62* | 0.58 | 0.55 | 0.71** |
| Salience | 0.48 | 0.51 | 0.49 | 0.52 |
| Probability | 0.54 | 0.61* | 0.57 | 0.53 |
| Identity | 0.45 | 0.47 | 0.46 | 0.44 |

`*` p < 0.05, `**` p < 0.01

**Note**: The exact layout will be finalized during implementation and this section updated accordingly.

---

### 9.5 Visualization Status Summary

| Figure | Script | Status |
|--------|--------|--------|
| Aggregated ANOVA | `plot_aggregated_anova.m` | ✓ Implemented |
| Aggregated Baseline | `plot_aggregated_baseline_comparison.m` | ✓ Implemented |
| Aggregated ROC | `plot_aggregated_roc_comparison.m` | ✓ Implemented |
| Aggregated Behavior | `plot_aggregated_behavior.m` | ✓ Implemented |
| Aggregated Decoding | `plot_aggregated_decoding.m` | △ Requires reorganization |
| ROC Scatter Plots | `plot_roc_scatter.m` | ○ Specified |
| SNc Subregion Comparison | `plot_snc_subregion.m` | ○ Specified |
| SC-SNc Interaction | `plot_sc_snc_interaction.m` | ○ Specified |
| Single-Neuron Summary PDFs | `generate_neuron_summary_pdf.m` | △ Requires revision |

**Status Legend:**
- ✓ Implemented: Functional and meeting current requirements
- △ Requires modification: Implemented but needs changes per this document
- ○ Specified: Defined in this document but not yet coded

**Note**: This table should be updated as implementation progresses.

---

## 10. Implementation Priorities

Implementation details (data structures, function signatures, code snippets) are maintained in a separate **Implementation Specification Document** to keep this rationale document focused on conceptual design.

### Phase 1: Core Infrastructure (Highest Priority)

| Task | Description |
|------|-------------|
| 1.1 | Add SNc subregion labeling to condition masks |
| 1.2 | Implement window-based ROC analysis |
| 1.3 | Create neuron metrics table structure for aggregation |
| 1.4 | Implement scatter plots with hover tooltips |

### Phase 2: Single-Neuron Summaries (High Priority)

| Task | Description |
|------|-------------|
| 2.1 | Redesign PDF to single-page layout (5 columns) |
| 2.2 | Add factor-specific PSTH panels |
| 2.3 | Integrate ROC summary table |

### Phase 3: Expanded Neuron Classification (Medium Priority)

| Task | Description |
|------|-------------|
| 3.1 | Implement neuron type classification (excited, suppressed, etc.) |
| 3.2 | Add suppressed neuron analysis for Hypothesis 4 |
| 3.3 | Update aggregation to include neuron types |

### Phase 4: SC-SNc Interaction (Medium Priority)

| Task | Description |
|------|-------------|
| 4.1 | Link simultaneous sessions via session_group_id |
| 4.2 | Implement factor-selectivity correlation analysis |
| 4.3 | Implement CCA analysis |

### Phase 5: Behavioral Analyses (Medium Priority)

| Task | Description |
|------|-------------|
| 5.1 | Add error rate analysis |
| 5.2 | Add peak velocity analysis (main sequence approach) |
| 5.3 | Add saccade direction as factor in behavioral ANOVA |

### Phase 6: Visualization Updates (Medium Priority)

| Task | Description |
|------|-------------|
| 6.1 | Reorganize decoding generalization plots (matrix layout) |
| 6.2 | Implement SNc subregion comparison figures |
| 6.3 | Implement SC-SNc interaction figures |

---

## Document History

| Date | Author | Description |
|------|--------|-------------|
| 2025-01-07 | Claude (with James Herman) | Initial creation based on codebase review and discussion |

---

## Notes

- Implementation details (code snippets, data structure specifications) are maintained in a separate Implementation Specification Document
- This document should be updated when analysis approaches change or when implementation reveals the need for design modifications
- The single-neuron PDF layout (Section 9.4) will be finalized during implementation
