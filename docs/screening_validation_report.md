# Screening Validation Report

**Generated:** 2026-02-04
**Analysis Script:** `code/analyze_screening_status.m`
**Data Source:** 11 sessions, 924 neurons

---

## Executive Summary

The neuron screening pipeline is fully implemented with all planned phases complete (Friedman caching, stable period detection, visualization). However, analysis reveals several issues requiring attention:

| Metric | Value | Assessment |
|--------|-------|------------|
| Overall inclusion rate | 37.2% | Moderate |
| SC inclusion rate | 38.8% | Moderate |
| SNc inclusion rate | 32.6% | Low |
| Session variability | 58.6% range | **High - Issue** |
| Neurons with unknown exclusion reason | 75 | **Bug** |
| Neurons near sparsity threshold | 148 SC | Sensitive |

**Primary Issues Identified:**
1. Bug: 75 excluded SC neurons have no exclusion reason set
2. High inter-session variability (12.5% to 71.1% inclusion)
3. Many neurons near sparsity threshold boundary

---

## 1. Current Screening Results

### Overall Statistics

```
Total neurons:           924
Total included:          344 (37.2%)

SC neurons:              694
SC included:             269 (38.8%)

SNc neurons:             230
SNc included:            75 (32.6%)
```

### Classification Breakdown

**SC Neurons:**
| Class | Count | Percentage |
|-------|-------|------------|
| classic | 138 | 19.9% |
| interneuron | 131 | 18.9% |
| excluded | 425 | 61.2% |

**SNc Neurons:**
| Class | Count | Percentage |
|-------|-------|------------|
| modulated | 75 | 32.6% |
| excluded | 155 | 67.4% |

### Per-Session Inclusion Rates

| Session | Area | Total | Included | Rate |
|---------|------|-------|----------|------|
| Feynman_08_12_2025_SNc | SNc | 50 | 27 | 54.0% |
| Feynman_08_12_2025_SC | SC | 56 | 17 | 30.4% |
| Feynman_08_13_2025_SNc | SNc | 38 | 27 | **71.1%** |
| Feynman_08_13_2025_SC | SC | 126 | 46 | 36.5% |
| Feynman_08_14_2025_SNc | SNc | 72 | 9 | **12.5%** |
| Feynman_08_14_2025_SC | SC | 99 | 21 | 21.2% |
| Feynman_08_15_2025_SNc | SNc | 70 | 12 | 17.1% |
| Feynman_08_15_2025_SC | SC | 111 | 46 | 41.4% |
| Newton_08_13_2025_SC | SC | 107 | 29 | 27.1% |
| Newton_08_14_2025_SC | SC | 103 | 48 | 46.6% |
| Newton_08_15_2025_SC | SC | 92 | 62 | **67.4%** |

**Observation:** Feynman_08_14_2025_SNc has anomalously low inclusion (12.5%), while Newton_08_15_2025_SC has high inclusion (67.4%). Investigate recording quality differences.

---

## 2. Exclusion Reasons Analysis

### Top Exclusion Reasons (n=580 excluded neurons)

| Count | Reason | Notes |
|-------|--------|-------|
| 207 | Low FR & Sparse | Combined failure |
| **75** | **(unknown)** | **BUG: Missing reason** |
| 68 | Low FR & No mod. & Sparse | Triple failure |
| 65 | No modulation | Modulation-only failure |
| 35 | No mod. | Abbreviated form |
| 19 | Low FR | FR-only failure |
| 11 | Sparse | Sparsity-only failure |

### Bug: Missing Exclusion Reasons

75 neurons (all SC) are classified as "excluded" but have empty `exclusion_reason` field. This occurs in `screen_sc_neurons.m` where the classification logic sets `neuron_class = 'excluded'` but doesn't populate the reason string.

**Location:** `screen_sc_neurons.m`, lines 614-617
**Fix Required:** Add exclusion reason assignment for SC neurons that fail classification.

---

## 3. Threshold Sensitivity

### Current Thresholds

| Parameter | SC | SNc |
|-----------|-----|-----|
| FR threshold | 3.0 sp/s | 1.0 sp/s |
| Sparsity threshold | 80% | 80% |

### SC Neuron Distribution

```
FR range:       0.01 - 430.00 sp/s (median: 4.02)
Sparsity range: 0.0% - 99.9% (median: 73.0%)

Marginal neurons (near threshold):
  - FR (within 20%):       67 neurons (9.7%)
  - Sparsity (within 10%): 148 neurons (21.3%)
```

### SNc Neuron Distribution

```
FR range:       0.02 - 71.19 sp/s (median: 8.38)
Sparsity range: 1.7% - 99.8% (median: 53.1%)

Marginal neurons (near threshold):
  - FR (within 20%):       4 neurons (1.7%)
  - Sparsity (within 10%): 26 neurons (11.3%)
```

### Yield Sensitivity Matrix (SC)

Shows % of neurons passing FR and sparsity thresholds (before modulation filter):

| FR \ Sparsity | 60% | 70% | 80% | 90% |
|---------------|-----|-----|-----|-----|
| 1.0 sp/s | 38.8% | 47.6% | 58.4% | 73.5% |
| 2.0 sp/s | 38.8% | 47.6% | 58.4% | 65.4% |
| 3.0 sp/s | 38.8% | 47.6% | **54.6%** | 56.3% |
| 4.0 sp/s | 38.8% | 47.1% | 49.7% | 50.3% |
| 5.0 sp/s | 38.8% | 45.1% | 45.8% | 46.0% |

**Current setting (3.0/80%)** yields 54.6% passing basic thresholds. Note:
- Relaxing sparsity from 80% to 90% only adds 1.7% yield
- Tightening FR from 3.0 to 4.0 loses 4.9% yield
- The 60% sparsity column shows a floor of ~39% regardless of FR

### Yield Sensitivity Matrix (SNc)

| FR \ Sparsity | 60% | 70% | 80% | 90% |
|---------------|-----|-----|-----|-----|
| 0.5 sp/s | 57.4% | 68.7% | 75.2% | 83.9% |
| 1.0 sp/s | 57.4% | 68.7% | **75.2%** | 83.9% |
| 2.0 sp/s | 57.4% | 68.7% | 75.2% | 78.3% |
| 3.0 sp/s | 57.4% | 68.7% | 73.0% | 73.0% |

**Observation:** SNc FR threshold has no effect between 0.5-1.5 sp/s because the median SNc FR (8.38 sp/s) is well above these values. Consider raising SNc FR threshold.

---

## 4. Stability Analysis

### Statistics

```
Neurons with stability info: 924/924 (100.0%)
Stability score:        median = 0.38, range = 0.00 - 1.00
Presence ratio:         median = 0.97, range = 0.00 - 1.00
Stable period fraction: median = 0.42, range = 0.01 - 1.00
```

### Interpretation

- **High presence ratio** (median 0.97): Most neurons fire in most time windows
- **Low stable period fraction** (median 0.42): Less than half of trials fall within detected stable periods
- **Moderate stability score** (median 0.38): Combined metric is lower due to stable period fraction

**Concern:** The stable period detection may be too strict, using only ~42% of trials for FR/sparsity calculation. This could artificially inflate FR estimates for neurons that were partially present.

---

## 5. Modulation Analysis

```
Total modulated: 651/924 (70.5%)
SC modulated:    522/694 (75.2%)
SNc modulated:   129/230 (56.1%)
```

### Observation

75.2% of SC neurons show task modulation, but only 38.8% are ultimately included. This suggests:
- FR and sparsity thresholds are the primary exclusion factors
- The spatial selectivity requirements for "classic" classification are restrictive

56.1% of SNc neurons show modulation, but only 32.6% pass all criteria, indicating FR/sparsity filtering removes ~23% of modulated neurons.

---

## 6. Identified Issues and Recommendations

### Issue 1: Missing Exclusion Reasons (Bug)

**Severity:** Medium
**Impact:** 75 neurons have unclear exclusion status
**Location:** `screen_sc_neurons.m`

**Recommendation:** Add exclusion reason assignment in the classification logic:
```matlab
if ~has_min_fr || is_sparse || ~is_task_modulated
    neuron_class{i_cluster} = 'excluded';
    selected_neurons(i_cluster) = false;
    % ADD: Build exclusion reason
    reasons = {};
    if ~has_min_fr; reasons{end+1} = 'Low FR'; end
    if is_sparse; reasons{end+1} = 'Sparse'; end
    if ~is_task_modulated; reasons{end+1} = 'No modulation'; end
    % Store in screening_info structure
end
```

### Issue 2: High Inter-Session Variability

**Severity:** High
**Impact:** Inclusion rates range from 12.5% to 71.1%

**Investigation Needed:**
1. Compare recording quality metrics (noise levels, spike sorting quality) between sessions
2. Check if low-yield sessions have different task conditions
3. Examine Feynman_08_14_2025_SNc specifically (12.5% yield)

**Possible Solutions:**
- Session-adaptive thresholds based on recording quality
- Normalize metrics within session before applying thresholds
- Manual review of outlier sessions

### Issue 3: Sparsity Threshold Sensitivity

**Severity:** Low-Medium
**Impact:** 148 SC neurons (21%) are within 10% of sparsity threshold

**Recommendation:** The current 80% threshold appears reasonable:
- Median SC sparsity is 73%, just below threshold
- Relaxing to 90% provides minimal benefit (+1.7%)
- Tightening to 70% would exclude 7% more neurons

No change recommended unless scientific rationale supports adjustment.

### Issue 4: SNc FR Threshold May Be Too Low

**Severity:** Low
**Impact:** Minimal - threshold has no effect at current level

**Observation:** SNc median FR is 8.38 sp/s, far above the 1.0 sp/s threshold. The yield is identical for thresholds 0.5-1.5 sp/s.

**Recommendation:** Consider raising SNc FR threshold to 2.0-3.0 sp/s to better discriminate quality neurons. However, this depends on biological expectations for SNc firing rates.

### Issue 5: Stable Period Detection Parameters

**Severity:** Medium
**Impact:** Only 42% of trials used for FR/sparsity calculation

**Investigation Needed:**
1. Review stability diagnostic plots to understand why stable periods are short
2. Check if neurons truly become unstable or if detection is too strict
3. Compare FR estimates with/without stable period restriction

**Possible Adjustments:**
- Increase `min_stable_trials` from 50 to require longer periods
- Decrease `cv_threshold` from 0.5 to be more lenient
- Add option to fall back to full session if stable period is < X% of trials

---

## 7. Recommended Action Plan

### Immediate (Bug Fixes)

1. [ ] Fix missing exclusion reasons in `screen_sc_neurons.m`
2. [ ] Re-run screening with `force_rerun.screening = true`
3. [ ] Regenerate screening metrics PDFs

### Short-Term (Validation)

4. [ ] Manually review 10 neurons from Feynman_08_14_2025_SNc (lowest yield session)
5. [ ] Manually review 10 neurons from Newton_08_15_2025_SC (highest yield session)
6. [ ] Compare classification against visual inspection
7. [ ] Review stability diagnostic plots for representative sessions

### Medium-Term (Parameter Tuning)

8. [ ] Investigate stable period detection parameters if Step 7 reveals issues
9. [ ] Consider raising SNc FR threshold if biologically appropriate
10. [ ] Document final threshold decisions in `analysis_rationale.md`

### Long-Term (Robustness)

11. [ ] Implement session quality metrics to explain variability
12. [ ] Consider session-normalized thresholds for future recordings
13. [ ] Create automated validation pipeline for new sessions

---

## 8. Files Modified/Created

| File | Status | Purpose |
|------|--------|---------|
| `docs/screening_finalization_plan.md` | Created | Implementation plan |
| `docs/screening_validation_report.md` | Created | This report |
| `code/analyze_screening_status.m` | Created | Analysis script |
| `data/screening_analysis_data.mat` | Created | Aggregated metrics |

---

## Appendix: Raw Session Data

```
Session                      | Area | n_total | n_included | pct_included | mean_fr | median_fr | pct_pass_fr | pct_pass_sparse | pct_modulated
Feynman_08_12_2025_SNc       | SNc  | 50      | 27         | 54.0%        | -       | -         | -           | -               | -
Feynman_08_12_2025_SC        | SC   | 56      | 17         | 30.4%        | -       | -         | -           | -               | -
Feynman_08_13_2025_SNc       | SNc  | 38      | 27         | 71.1%        | -       | -         | -           | -               | -
Feynman_08_13_2025_SC        | SC   | 126     | 46         | 36.5%        | -       | -         | -           | -               | -
Feynman_08_14_2025_SNc       | SNc  | 72      | 9          | 12.5%        | -       | -         | -           | -               | -
Feynman_08_14_2025_SC        | SC   | 99      | 21         | 21.2%        | -       | -         | -           | -               | -
Feynman_08_15_2025_SNc       | SNc  | 70      | 12         | 17.1%        | -       | -         | -           | -               | -
Feynman_08_15_2025_SC        | SC   | 111     | 46         | 41.4%        | -       | -         | -           | -               | -
Newton_08_13_2025_SC         | SC   | 107     | 29         | 27.1%        | -       | -         | -           | -               | -
Newton_08_14_2025_SC         | SC   | 103     | 48         | 46.6%        | -       | -         | -           | -               | -
Newton_08_15_2025_SC         | SC   | 92      | 62         | 67.4%        | -       | -         | -           | -               | -
```

(Full per-session metrics available in `data/screening_analysis_data.mat`)
