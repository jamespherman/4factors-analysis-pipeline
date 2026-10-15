# SC Neuron Population Report

_Generated: 2026-02-01 10:56:24_

## Summary

- **Total SC neurons recorded**: 694
- **Total selected (classic + interneuron)**: 269
- **Classic SC neurons**: 138
- **Putative interneurons**: 131
- **Excluded**: 425
- **Total sessions**: 7

## Neuron Counts by Subject

| Subject | Classic | Interneuron | Total Selected | Excluded | Total Recorded | Sessions |
|---------|---------|-------------|----------------|----------|----------------|----------|
| Feynman | 57 | 73 | 130 | 262 | 392 | 4 |
| Newton | 81 | 58 | 139 | 163 | 302 | 3 |
| **Total** | **138** | **131** | **269** | **425** | **694** | **7** |

## Neuron Counts by Session

| Session | Subject | Date | Classic | Interneuron | Total Selected | Excluded | Total Recorded |
|---------|---------|------|---------|-------------|----------------|----------|----------------|
| Feynman_08_12_2025_SC | Feynman | 08/12/2025 | 12 | 5 | 17 | 39 | 56 |
| Feynman_08_13_2025_SC | Feynman | 08/13/2025 | 8 | 38 | 46 | 80 | 126 |
| Feynman_08_14_2025_SC | Feynman | 08/14/2025 | 14 | 7 | 21 | 78 | 99 |
| Feynman_08_15_2025_SC | Feynman | 08/15/2025 | 23 | 23 | 46 | 65 | 111 |
| Newton_08_13_2025_SC | Newton | 08/13/2025 | 15 | 14 | 29 | 78 | 107 |
| Newton_08_14_2025_SC | Newton | 08/14/2025 | 22 | 26 | 48 | 55 | 103 |
| Newton_08_15_2025_SC | Newton | 08/15/2025 | 44 | 18 | 62 | 30 | 92 |

## Target Onset Responses by Location (Visual Epoch)

Counts of neurons with significant activation or suppression at target onset, by target location.

| Location | Angle | Classic Activated | Classic Suppressed | Interneuron Activated | Interneuron Suppressed |
|----------|-------|-------------------|--------------------|-----------------------|------------------------|
| 1 | 45° | 2 (4.7%) | 2 (4.7%) | 15 (22.7%) | 7 (10.6%) |
| 2 | 135° | 25 (37.9%) | 0 (0.0%) | 36 (40.4%) | 8 (9.0%) |
| 3 | 225° | 12 (27.9%) | 0 (0.0%) | 29 (43.9%) | 7 (10.6%) |
| 4 | 315° | 0 (0.0%) | 2 (4.7%) | 25 (37.9%) | 10 (15.2%) |

_Note: Each neuron is counted once per location where it was tested. Percentages are relative to the number of neurons tested at each location._

## 4-Factors Modulation (Visual Epoch)

Proportion of neurons with significant (p < 0.05) selectivity for each factor during the visual response epoch.

| Factor | Classic Sig | Classic Total | Classic % | Interneuron Sig | Interneuron Total | Interneuron % |
|--------|-------------|---------------|-----------|-----------------|-------------------|---------------|
| Reward | 5 | 119 | 4.2% | 3 | 114 | 2.6% |
| Salience | 10 | 119 | 8.4% | 6 | 113 | 5.3% |
| Identity | 1 | 119 | 0.8% | 2 | 113 | 1.8% |
| Probability | 1 | 117 | 0.9% | 1 | 114 | 0.9% |

## Methods

### Neuron Classification

**Classic SC neurons** meet all criteria:
- Significant increase in visual/delay/saccade epoch vs baseline for ≥1 contralateral location
- Mean FR > 3.0 sp/s
- NOT significantly activated at opposing locations (180° apart)
- Not sparse (≤80% of 100ms bins empty)

**Putative interneurons** meet:
- Mean FR > 3.0 sp/s
- Task-modulated (activation or suppression somewhere)
- Either opposing activation OR contralateral suppression

**Excluded** neurons fail to meet basic criteria (low FR, sparse, or not task-modulated).

### Statistical Tests

- Target onset responses: Friedman test with post-hoc Wilcoxon signed-rank tests (Bonferroni corrected)
- 4-Factors modulation: Window-based ROC analysis comparing high vs low conditions
