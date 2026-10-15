# Screening Finalization Plan

**Created:** 2026-02-04
**Purpose:** Finalize neuron screening criteria for SC and SNc data
**Status:** Active implementation plan

---

## Executive Summary

Based on comprehensive research of the codebase, the screening infrastructure is **fully implemented**:
- Phase 1 (Friedman caching): Complete
- Phase 2 (Stable period detection): Complete
- Phase 3 (Visualization): Complete and integrated

The task is now **validation and tuning**, not implementation. This plan outlines the steps to:
1. Validate existing outputs
2. Assess screening effectiveness
3. Identify specific issues
4. Iterate on parameters
5. Finalize and document criteria

---

## Current Implementation Status

### Files Implemented and Integrated

| Component | File | Status |
|-----------|------|--------|
| Screening orchestrator | `apply_neuron_screening.m` | Complete |
| Stable period detection | `detect_stable_trial_period.m` | Complete |
| SC classification | `screen_sc_neurons.m` | Complete |
| SNc classification | `screen_snc_neurons.m` | Complete |
| Session visualization | `plot_screening_metrics_session.m` | Complete |
| Aggregate visualization | `plot_screening_metrics_aggregate.m` | Complete |
| Data aggregation | `aggregate_screening_data.m` | Complete |
| Pipeline integration | `run_4factors_analysis.m` (lines 478, 963-964) | Complete |

### Current Thresholds

**SC Neurons (`screen_sc_neurons.m`):**
- `MIN_FR_THRESHOLD = 3.0` sp/s
- `MAX_SPARSE_FRACTION = 0.80` (80% empty bins)
- Classification: classic, interneuron, excluded

**SNc Neurons (`screen_snc_neurons.m`):**
- `fr_threshold = 1.0` sp/s (from `define_task_conditions.m`)
- `max_proportion_empty = 0.80` (80% empty bins)
- Classification: modulated, excluded

**Stable Period Detection (`detect_stable_trial_period.m`):**
- `window_size = 20` trials
- `step_size = 5` trials
- `cv_threshold = 0.5`
- `min_fr_threshold = 0.1` sp/s
- `min_stable_trials = 50`

---

## Step-by-Step Execution Plan

### Step 1: Generate Fresh Visualizations

**Objective:** Ensure all screening metrics PDFs are current and complete.

**Actions:**
1. Set `force_rerun.screening = true` in `run_4factors_analysis.m`
2. Ensure `generate_screening_metrics_plots = true`
3. Run the pipeline: `matlab -batch "run('code/run_4factors_analysis.m')"`
4. Verify outputs in `figures/`:
   - `{session_id}_screening_metrics.pdf` for each session
   - `aggregate_screening_metrics.pdf`

**Success Criteria:**
- All session PDFs generated without errors
- Aggregate PDF contains all 6 planned figures

**Estimated Duration:** Pipeline runtime (~30-60 min depending on caching)

---

### Step 2: Review Session-Level Visualizations

**Objective:** Identify per-session screening anomalies.

**For each session, examine `{session_id}_screening_metrics.pdf`:**

**Page 1 - FR vs Sparsity Scatter:**
- [ ] Are threshold lines positioned correctly?
- [ ] How many neurons fall near threshold boundaries?
- [ ] Are there obvious outliers?
- [ ] Is the included/excluded coloring sensible?

**Page 2 - Stability Metrics:**
- [ ] What fraction of neurons have stable periods < full session?
- [ ] Is stable FR vs total FR correlation reasonable (near diagonal)?
- [ ] Are stability scores distributed sensibly?

**Page 3 - Modulation Summary:**
- [ ] Is modulation detection working (non-zero counts by event)?
- [ ] Are exclusion reasons distributed as expected?

**Page 4 - Exclusion Breakdown:**
- [ ] What's the primary exclusion reason?
- [ ] Are there unexpected patterns?

**Output:** Notes file documenting issues per session.

---

### Step 3: Review Aggregate Visualizations

**Objective:** Identify population-level patterns and cross-session consistency.

**Examine `aggregate_screening_metrics.pdf`:**

**Figure 2.1 - Screening Yield Summary:**
- [ ] Are yields consistent across sessions?
- [ ] Any sessions with unusually high/low inclusion rates?

**Figure 2.2 - Population Metric Distributions:**
- [ ] Do SC and SNc distributions differ as expected?
- [ ] Are thresholds positioned at sensible points in distributions?

**Figure 2.3 - Aggregate Scatter Plot:**
- [ ] Is there clear separation between included/excluded?
- [ ] Are there clusters of near-threshold neurons?

**Figure 2.5 - Threshold Sensitivity Analysis:**
- [ ] How steep is the yield curve at current thresholds?
- [ ] Would small threshold changes dramatically alter yield?

**Output:** Summary document with specific concerns.

---

### Step 4: Threshold Sensitivity Analysis

**Objective:** Quantify impact of threshold changes on screening yield.

**Actions:**
1. Run `analyze_threshold_sensitivity.m` (if exists) or create ad-hoc analysis
2. Test FR threshold range: [1.0, 2.0, 3.0, 4.0, 5.0] sp/s for SC
3. Test sparsity threshold range: [0.60, 0.70, 0.80, 0.90]
4. Document yield at each combination

**Key Questions:**
- What threshold combination gives ~50% inclusion rate?
- What threshold combination matches manual inspection expectations?
- Are any neurons "marginal" (within 10% of thresholds)?

**Output:** Threshold sensitivity table with recommendations.

---

### Step 5: Manual Neuron Inspection

**Objective:** Validate automated classification against expert judgment.

**Sample Selection:**
1. Select 5 neurons classified as "classic SC" - verify spatial selectivity
2. Select 5 neurons classified as "interneuron" - verify opposing activation
3. Select 5 neurons classified as "excluded" - verify they should be excluded
4. Select 5 "marginal" neurons (near thresholds) - determine correct classification

**For each neuron, examine:**
- Neuron diagnostic PDF (if exists)
- Spatial tuning plots
- PSTH by condition
- Stability raster

**Output:** Validation spreadsheet with manual vs automated classification agreement.

---

### Step 6: Parameter Adjustment (if needed)

**Objective:** Update thresholds based on validation findings.

**If adjustments needed:**
1. Update thresholds in `define_task_conditions.m` (neuron_inclusion struct)
2. For SC-specific: update constants in `screen_sc_neurons.m`
3. Re-run screening with `force_rerun.screening = true`
4. Re-generate visualizations
5. Repeat Steps 2-5 until satisfied

**Documentation:**
- Record rationale for each threshold change
- Save before/after comparison metrics

---

### Step 7: Finalize and Document

**Objective:** Lock down final criteria and document decisions.

**Deliverables:**
1. Update `define_task_conditions.m` with final thresholds
2. Update `docs/design/analysis_rationale.md` with screening justification
3. Create summary table of final screening statistics:
   - Total neurons per session
   - Neurons passing each criterion
   - Final yield by brain area
4. Commit changes with descriptive message

---

## Risk Assessment

### Known Issues to Watch For

1. **Stable period detection may be too strict**
   - If many neurons have `is_stable = false`, parameters may need loosening
   - Check `presence_ratio` distribution

2. **SC opposing activation criterion may be too sensitive**
   - Some classic neurons may be misclassified as interneurons
   - Manual overrides may be needed (as already done for Newton_08_14_2025_SC cluster 6)

3. **SNc modulation detection may miss slow responses**
   - 3-consecutive-bins requirement may be too strict for slow SNc dynamics
   - Consider reviewing modulation windows

4. **Cross-session variability**
   - Different recording quality may require session-specific thresholds
   - Currently uses same thresholds for all sessions

---

## Execution Commands

```matlab
% Step 1: Run full pipeline with fresh screening
cd('D:\OneDrive - University of Pittsburgh\Code\4factors-analysis-pipeline')
run('code/run_4factors_analysis.m')

% Alternative: Run just screening on one session for testing
session_data = load('path/to/session_data.mat');
[~, analysis_plan] = define_task_conditions();
[selected_neurons, screening_info, neuron_class, session_data] = ...
    apply_neuron_screening(session_data, analysis_plan, project_root);

% Generate just the aggregate visualization
aggregate_data = aggregate_screening_data(manifest, project_root);
plot_screening_metrics_aggregate(aggregate_data, fullfile(project_root, 'figures'));
```

---

## Success Criteria

Screening will be considered "finalized" when:

1. [ ] All session visualizations reviewed without critical issues
2. [ ] Aggregate statistics are documented
3. [ ] Manual inspection validates automated classification (>90% agreement)
4. [ ] Threshold sensitivity is understood and documented
5. [ ] Final thresholds are recorded in `define_task_conditions.m`
6. [ ] Rationale is documented in `analysis_rationale.md`
7. [ ] No known bugs or anomalies in screening logic

---

## Autonomous Execution Notes

**Steps that can be fully automated:**
- Step 1: Pipeline execution
- Partial Step 4: Threshold sensitivity computation

**Steps requiring human judgment:**
- Steps 2-3: Visual review of PDFs
- Step 5: Manual neuron inspection
- Step 6: Decision on parameter changes
- Step 7: Final approval

**Recommended approach:**
- Claude can execute Step 1 and prepare summaries
- Claude can identify potential issues in visualizations
- Human reviews and approves changes
- Claude implements approved parameter updates

---

## Appendix: Key File Locations

```
code/
├── apply_neuron_screening.m      # Main screening orchestrator
├── screen_sc_neurons.m           # SC classification logic
├── screen_snc_neurons.m          # SNc classification logic
├── detect_stable_trial_period.m  # Stability detection
├── compute_task_modulation.m     # arrayROC-based modulation
├── plot_screening_metrics_session.m   # Session-level plots
├── plot_screening_metrics_aggregate.m # Population-level plots
├── aggregate_screening_data.m    # Data collection utility
└── utils/
    └── define_task_conditions.m  # Threshold definitions

figures/
├── {session_id}_screening_metrics.pdf    # Per-session outputs
├── aggregate_screening_metrics.pdf       # Population summary
└── {session_id}_snc_screening.pdf        # SNc-specific diagnostics

docs/
├── screening_improvement_plan.md         # Original implementation plan
├── screening_visualizations_plan.md      # Visualization specifications
└── screening_finalization_plan.md        # THIS DOCUMENT
```
