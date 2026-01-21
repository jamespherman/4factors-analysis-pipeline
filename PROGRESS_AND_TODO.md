# Project Progress and To-Do List

This document tracks the development of the 4factors analysis pipeline.

**Last Updated:** 2026-01-21 (Task 8.3 SNc classification with arrayROC-based modulation detection)

---

## Status Legend

- [x] Complete
- [~] Partially complete / Needs revision
- [ ] Not started

---

## Part 0: Project Setup and Repository Initialization

- [x] Create GitHub repository with standard MATLAB .gitignore
- [x] Create initial directory structure: /code, /config, /docs, /output
- [x] Add helper functions to `/code/utils/`
- [x] Add README.md and PROGRESS_AND_TODO.md
- [x] Create `config/session_manifest.csv` with status columns

---

## Part 1: Core Pipeline Implementation

- [x] **Task 1.1: Master Script (`run_4factors_analysis.m`)**
  - Reads session manifest, loads analysis plan, iterates sessions, saves results

- [x] **Task 1.2: Neuron Screening**
  - `screen_da_neurons.m`: Selects putative DA neurons (firing rate + waveform)
  - `screen_sc_neurons.m`: Identifies task-modulated SC neurons

- [x] **Task 1.3: Data Preparation (`prepare_core_data.m`)**
  - Aligns and bins spike data for all events/conditions

- [x] **Task 1.4: Analysis Plan (`define_task_conditions.m`)**
  - Returns `analysis_plan` struct (no args) or condition masks (with session_data)

---

## Part 2: Single-Neuron Analysis

- [x] **Task 2.1: Neuron Diagnostic PDFs (`generate_neuron_summary_pdf.m`)**
  - Generates diagnostic plots for each selected neuron

- [x] **Task 2.2: Baseline Comparison (`analyze_baseline_comparison.m`)**
  - ROC analysis comparing post-event to pre-event baseline

- [x] **Task 2.3: Time-Resolved ROC (`analyze_roc_comparison.m`)**
  - Sliding window ROC between experimental conditions

- [x] **Task 2.4: N-way ANOVA (`analyze_anova.m`)**
  - Main effects and interactions of task factors

---

## Part 3: Aggregation and Population Plotting

- [x] **Task 3.1: Aggregate Results (`aggregate_analysis_results.m`)**
  - Concatenates single-session results into aggregated file

- [x] **Task 3.2: Plot Aggregated Results**
  - `plot_aggregated_baseline_comparison.m`
  - `plot_aggregated_roc_comparison.m`
  - `plot_aggregated_anova.m`
  - `plot_aggregated_behavior.m`

---

## Part 4: Two-Model Refactoring (Confounded Factors)

- [x] **Task 4.1: Split Analysis Plans**
  - Updated `define_task_conditions.m` for Image vs Bullseye models

- [x] **Task 4.2: Update Aggregation**
  - `aggregate_analysis_results.m` handles separate model structures

- [x] **Task 4.3: Update Plotting**
  - Visualization for two-model results

---

## Part 5: Population Decoding

- [x] **Task 5.1: SVM Classification**
  - Time-resolved and windowed decoding of task variables

- [~] **Task 5.2: Decoding Visualization (`plot_aggregated_decoding.m`)**
  - Current: Flat grid layout
  - **Needs revision**: Reorganize into matrix layout (train × test factors)

---

## Part 6: SNc Subregion Analysis (NEW)

- [x] **Task 6.1: Subregion Labeling**
  - Parse grid coordinates from manifest
  - Add `snc_subregion` field to session metadata
  - Create `is_rvmSNc` and `is_cdlSNc` condition masks
  - Rule: Y ∈ {3,4} → rvmSNc; Y ∈ {1,2} → cdlSNc (SNc sessions only)

- [x] **Task 6.2: Subregion Comparison Figures**
  - Grouped bar plot: mean ROC by factor × subregion
  - Violin plots: ROC distributions by subregion
  - Statistical annotations (Wilcoxon rank-sum)
  - Implementation: `plot_snc_subregion_comparison.m`

---

## Part 7: Window-Based ROC and Neuron Metrics (NEW)

- [x] **Task 7.1: Window-Based ROC (`analyze_window_roc.m`)**
  - Fixed-window ROC AUC (not max-in-window)
  - Area-specific windows: SC visual [50, 250]ms; SNc visual [100, 300]ms
  - Output: Single summary metric per neuron/factor/epoch

- [x] **Task 7.2: Neuron Metrics Table**
  - Create `neuron_metrics_table` in aggregation
  - Fields: session_id, cluster_id, brain_area, snc_subregion, neuron_type, ROC values per epoch/factor

- [x] **Task 7.3: ROC Scatter Plots (`plot_roc_scatter.m`)**
  - Pairwise factor comparisons (6 panels for 4 factors)
  - Hover tooltips for neuron identification
  - CSV backup for each plot

---

## Part 8: Expanded Neuron Classification (NEW)

- [x] **Task 8.1: SC Neuron Classification Categories**
  - Expand beyond binary `selected_neurons`
  - Categories: 'classic', 'interneuron', 'excluded'
  - Integrated: `apply_neuron_screening.m` calls `screen_sc_neurons.m` for SC sessions
  - `neuron_class` stored in `screening_info` and propagated to `neuron_type` in aggregation

- [ ] **Task 8.2: Putative SC Interneuron Neuron Analysis**
  - Test Hypothesis 4: Putative SC interneurons show factor-non-selective responses
  - Compute factor-specific ROC for Putative SC interneuron neurons
  - Compare to chance (ROC ≈ 0.5)

- [x] **Task 8.3: SNc Neuron Classification Categories**
  - Categories: 'modulated', 'excluded'
  - Implementation: `screen_snc_neurons.m` called from `apply_neuron_screening.m`
  - Criteria for 'modulated' (ALL must be true):
    - (a) Significant change in activity (increase OR decrease) for at least 3 consecutive
          samples using arrayROC bootstrap CI (95% excludes 0.5). Checked at: fixOn,
          targetOn, fixOff, saccadeOnset, reward via `compute_task_modulation.m`
    - (b) Mean FR > 1 sp/s (computed over active period)
    - (c) Not sparse (<=70% of 100ms bins empty across active period)
  - Generates diagnostic figure: `{session_id}_snc_screening.pdf`

- [x] **Task 8.4: Update Aggregation**
  - Include `neuron_type` categorical variable in `neuron_metrics_table`
  - Added screening metrics: `screening_mean_fr`, `screening_sparsity`,
    `passes_fr_threshold`, `passes_sparsity_threshold`, `is_task_modulated`
  - Filter analyses by neuron type (enabled by `neuron_type` column)

---

## Part 9: Enhanced Neuron Summary PDFs (NEW)

- [x] **Task 9.1: Single-Page Layout Redesign**
  - 5-column layout: 2 spatial (jph task) + 3 factor (targetOn, saccadeOnset, reward)
  - Include: waveform, ISI, spatial tuning, factor PSTHs, ROC summary table
  - Scope: Generate for modulated + suppressed neurons (not unmodulated)

- [x] **Task 9.2: ROC Summary Table Integration**
  - Embed factor × epoch ROC table with significance markers
  - Enable verification against scatter plot positions
  - Displays window-based ROC AUC values with * for p < 0.05

---

## Part 10: SC-SNc Interaction Analysis (NEW)

- [ ] **Task 10.1: Link Simultaneous Sessions**
  - Use `session_group_id` from manifest to identify paired SC+SNc recordings
  - Available: 4 sessions (2 rvmSNc, 2 cdlSNc)

- [ ] **Task 10.2: Factor Selectivity Correlation**
  - Correlate SC population decoding with SNc subregion-specific ROC
  - Separate analysis for stimulus-driven vs goal-directed factors

- [ ] **Task 10.3: CCA Analysis (Optional)**
  - Canonical Correlation Analysis between SC and SNc populations
  - Test: Do CCA weights correlate with factor selectivity by subregion?

---

## Part 11: Behavioral Analysis Enhancements (NEW)

- [ ] **Task 11.1: Error Rate Analysis**
  - Proportion of non-rewarded trials by condition
  - Logistic regression or chi-square test

- [ ] **Task 11.2: Saccade Direction Effects**
  - Add direction (left/right) as factor in behavioral ANOVA
  - Test for direction × factor interactions

- [ ] **Task 11.3: Main Sequence Analysis**
  - Analyze peak velocity directly (amplitude constant across conditions)
  - Optional: include amplitude as covariate for landing variability

---

## Part 12: Independence Testing (NEW)

- [ ] **Task 12.1: Chi-Square Independence Test**
  - Test whether dual-factor selectivity differs from independence
  - Expected under independence: P(A∩B) = P(A) × P(B)
  - Apply to all factor pairs

---

## Part 13: Documentation (NEW)

- [x] **Task 13.1: Analysis Rationale Document**
  - Scientific hypotheses, analysis logic, interpretation
  - Location: `docs/design/analysis_rationale.md`

- [~] **Task 13.2: Implementation Specification Document**
  - Data structures, function signatures, code patterns
  - Location: `docs/design/implementation_spec.md`
  - Status: Starter created, to be populated during implementation

- [x] **Task 13.3: Update AGENTS.md**
  - Revised for Claude Code workflow
  - Added references to design documents

---

## Deferred / Removed

- ~~**Dimensionality Reduction (dPCA)**~~
  - Removed from scope per design discussion
  - May revisit if core analyses are complete

---

## Implementation Priority Order

Based on analysis rationale document:

### Phase 1: Core Infrastructure (Current Focus)
1. Task 6.1: SNc subregion labeling
2. Task 7.1: Window-based ROC
3. Task 7.2: Neuron metrics table
4. Task 7.3: ROC scatter plots with tooltips

### Phase 2: Enhanced Summaries
1. Task 9.1: Single-page PDF layout
2. Task 9.2: ROC summary table integration

### Phase 3: Expanded Classification
1. Task 8.1: SC Classification categories (classic vs interneuron)
2. Task 8.2: Putative interneuron analysis
3. Task 8.3: SNc Classification categories
4. ~~Task 8.4: Update aggregation~~ (COMPLETE)

### Phase 4: SC-SNc Interaction
1. Task 10.1: Link simultaneous sessions
2. Task 10.2: Factor selectivity correlation
3. Task 10.3: CCA analysis (optional)

### Phase 5: Behavioral & Statistical
1. Task 11.1-11.3: Behavioral enhancements
2. Task 12.1: Chi-square independence test

### Phase 6: Visualization Updates
1. Task 5.2: Decoding plot reorganization
2. Task 6.2: SNc subregion comparison figures
