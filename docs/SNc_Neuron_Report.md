# SNc Neuron Population Report

_Generated: 2026-02-01 13:18:53_

## Summary

- **Total SNc neurons recorded**: 230
- **Selected (modulated)**: 75 (32.6%)
- **Excluded**: 155 (67.4%)
- **Total sessions**: 4

## Neuron Counts by Subject

| Subject | Modulated | Excluded | Total | Sessions |
|---------|-----------|----------|-------|----------|
| Feynman | 75 | 155 | 230 | 4 |

## Neuron Counts by Session

| Session | Subject | Date | Modulated | Excluded | Total |
|---------|---------|------|-----------|----------|-------|
| Feynman_08_12_2025_SNc | Feynman | 08/12/2025 | 27 | 23 | 50 |
| Feynman_08_13_2025_SNc | Feynman | 08/13/2025 | 27 | 11 | 38 |
| Feynman_08_14_2025_SNc | Feynman | 08/14/2025 | 9 | 63 | 72 |
| Feynman_08_15_2025_SNc | Feynman | 08/15/2025 | 12 | 58 | 70 |

## Baseline Firing Rate Distribution (Selected Neurons)

Distribution of baseline firing rates among the 75 selected (modulated) SNc neurons.

| Baseline FR Range | Count | Percentage |
|-------------------|-------|------------|
| < 5 sp/s | 20 | 26.7% |
| 5-10 sp/s | 20 | 26.7% |
| 10-20 sp/s | 20 | 26.7% |
| >= 20 sp/s | 15 | 20.0% |
| **< 20 sp/s (Total)** | **60** | **80.0%** |

## Baseline Firing Rate by Session (Selected Neurons)

| Session | Modulated | < 20 sp/s | >= 20 sp/s | % Low FR |
|---------|-----------|-----------|------------|----------|
| Feynman_08_12_2025_SNc | 27 | 20 | 7 | 74.1% |
| Feynman_08_13_2025_SNc | 27 | 23 | 4 | 85.2% |
| Feynman_08_14_2025_SNc | 9 | 7 | 2 | 77.8% |
| Feynman_08_15_2025_SNc | 12 | 10 | 2 | 83.3% |

## Baseline FR Summary Statistics (Selected Neurons)

| Statistic | Value |
|-----------|-------|
| Mean | 11.55 sp/s |
| Median | 8.94 sp/s |
| Std Dev | 9.86 sp/s |
| Min | 0.00 sp/s |
| Max | 45.47 sp/s |

## Methods

### SNc Neuron Classification

**Modulated SNc neurons** meet all criteria:
- Significant change in activity (increase or decrease) for at least 3 consecutive samples
- Modulation detected at any of: fixation onset, target onset, fixation offset, peri-saccade, or reward
- Mean FR > 1 sp/s (computed over active period)
- Not sparse (<=80% of 100ms bins empty)

**Excluded** neurons fail to meet one or more of the above criteria.

### Baseline Firing Rate

Baseline FR is computed during the pre-target baseline epoch (-75 to +25 ms relative to target onset).
The 20 sp/s threshold is used to identify putative dopamine neurons, as dopamine neurons typically
have low baseline firing rates (<10 sp/s), while GABAergic interneurons often have higher rates.
