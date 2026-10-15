# Screening Visualizations Plan

## Overview

This plan outlines visualizations to understand the conservativeness/liberalness of neuron screening criteria. Visualizations are organized into:
1. **Session-level** - Per-session diagnostic plots showing individual neuron distributions
2. **Across-sessions (aggregate)** - Population-level summaries combining all sessions

## Goals

- Understand the distribution of key screening metrics (FR, sparsity, stability)
- Visualize where thresholds fall relative to the population
- Compare included vs. excluded neurons
- Identify sessions with unusual screening characteristics
- Track screening yield across recording sessions

---

## Part 1: Session-Level Visualizations

### Figure 1.1: Screening Metrics Scatter Plot (Primary Diagnostic)

**Purpose:** Show the 2D relationship between firing rate and sparsity with threshold boundaries.

**Layout:** Single figure with threshold lines and classification coloring

**Content:**
- X-axis: Mean firing rate (log scale, sp/s)
- Y-axis: Proportion of empty bins (0-1)
- Vertical line: FR threshold (3.5 sp/s for SC, 1.0 sp/s for SNc)
- Horizontal line: Sparsity threshold (0.70)
- Point color: Classification (included/excluded, or by exclusion reason)
- Point size: Stability score (larger = more stable)
- Annotations: Count in each quadrant (pass both, fail FR only, fail sparsity only, fail both)

**Key insight:** Shows how many neurons are near the thresholds vs. clearly passing/failing.

---

### Figure 1.2: Stability Metrics Panel

**Purpose:** Show the distribution of stability-related metrics.

**Layout:** 2x2 subplot panel

**Subplots:**
1. **Histogram: Stable period duration** (as % of total trials)
   - X: Proportion of session covered by stable period
   - Vertical line: Current minimum requirement

2. **Histogram: Presence ratio**
   - X: Fraction of windows with meaningful firing
   - Compare included vs. excluded neurons (overlapping histograms)

3. **Histogram: Stability score**
   - X: Combined stability metric (0-1)
   - Color by included/excluded

4. **Scatter: Stable FR vs. Total FR**
   - X: Total session FR
   - Y: Stable period FR
   - Diagonal line: y=x (no change)
   - Shows neurons where FR estimate changes significantly with stability correction

---

### Figure 1.3: Modulation Summary

**Purpose:** Show task modulation characteristics.

**Layout:** 2x2 subplot panel

**Subplots:**
1. **Pie chart: Modulation status**
   - Segments: Modulated (included), Modulated (excluded for other reasons), Not modulated

2. **Bar chart: Modulation by event type**
   - X: Event types (targetOn, saccadeOnset, etc.)
   - Y: Number of neurons showing modulation at each event
   - Stacked by direction (increase/decrease)

3. **Histogram: Number of modulated locations per neuron**
   - X: Count of locations showing modulation (0 to max)
   - Color by included/excluded

4. **Text summary panel**
   - Total neurons
   - Pass FR: N (X%)
   - Pass sparsity: N (X%)
   - Pass modulation: N (X%)
   - Final included: N (X%)

---

### Figure 1.4: Exclusion Breakdown (Sankey-style or Stacked Bar)

**Purpose:** Show the flow from total neurons to final selection through each filter.

**Layout:** Single horizontal stacked bar or waterfall chart

**Content:**
- Starting: Total neurons
- After FR filter: N remaining
- After sparsity filter: N remaining
- After modulation filter: N remaining
- Final: Selected neurons

Alternative: Venn-style diagram showing overlap of exclusion reasons.

---

## Part 2: Across-Sessions Visualizations

### Figure 2.1: Screening Yield Summary

**Purpose:** Compare screening outcomes across all sessions.

**Layout:** Grouped bar chart or heatmap

**Content:**
- X-axis: Session ID (grouped by subject, then by brain area)
- Y-axis: Percentage of neurons
- Bars/cells:
  - % passing FR threshold
  - % passing sparsity threshold
  - % showing modulation
  - % final included
- Color intensity: Yield percentage

**Key insight:** Identify sessions with unusually high/low yields.

---

### Figure 2.2: Population Metric Distributions

**Purpose:** Show overall distributions combining all neurons across sessions.

**Layout:** 3x2 subplot panel

**Subplots:**
1. **Histogram: Mean firing rate (all neurons)**
   - Log-scale X-axis
   - Vertical lines: SC threshold (3.5), SNc threshold (1.0)
   - Separate colors for SC vs. SNc neurons

2. **Histogram: Sparsity (all neurons)**
   - Vertical line: Threshold (0.70)
   - Separate colors for SC vs. SNc

3. **Histogram: Stability score (all neurons)**
   - Compare included vs. excluded

4. **Histogram: Stable period duration (% of trials)**
   - Compare included vs. excluded

5. **Histogram: Presence ratio**
   - Compare included vs. excluded

6. **Box plots: FR by session**
   - X: Session (ordered by date)
   - Y: Mean FR
   - Separate box for included vs. all neurons

---

### Figure 2.3: Aggregate Scatter Plot (Grand Summary)

**Purpose:** Show all neurons from all sessions in a single FR vs. sparsity space.

**Layout:** Single large scatter plot

**Content:**
- X: Mean firing rate (log scale)
- Y: Proportion empty bins
- Color: Included (green) vs. Excluded (red/gray)
- Shape: Brain area (circle=SC, triangle=SNc)
- Transparency: Alpha=0.5 to show density
- Marginal histograms: On top and right showing distributions
- Threshold lines: FR and sparsity thresholds
- Density contours: Optional, to show population density

---

### Figure 2.4: Session Comparison Panel

**Purpose:** Allow direct comparison of screening metrics between sessions.

**Layout:** Multi-panel figure with one row per session

**Content per row:**
- Mini FR histogram (normalized)
- Mini sparsity histogram (normalized)
- Yield bar (% included)
- N neurons text

**Key insight:** Quickly identify outlier sessions.

---

### Figure 2.5: Threshold Sensitivity Analysis

**Purpose:** Show how yield would change with different thresholds.

**Layout:** 2x2 subplot panel

**Subplots:**
1. **Line plot: Yield vs. FR threshold**
   - X: FR threshold (0.5 to 5.0 sp/s)
   - Y: % neurons included
   - Separate lines for SC and SNc
   - Vertical line: Current threshold

2. **Line plot: Yield vs. sparsity threshold**
   - X: Sparsity threshold (0.5 to 0.9)
   - Y: % neurons included
   - Vertical line: Current threshold

3. **Heatmap: Yield vs. both thresholds**
   - X: FR threshold
   - Y: Sparsity threshold
   - Color: % included
   - Star marker: Current thresholds

4. **Text summary**
   - Current thresholds and yields
   - Marginal neurons (within 10% of thresholds)

---

### Figure 2.6: Temporal Trends

**Purpose:** Check for drift in recording quality or screening outcomes over time.

**Layout:** 2x1 subplot panel

**Subplots:**
1. **Line plot: Metrics over sessions**
   - X: Session date
   - Y: Mean values (separate lines for mean FR, mean sparsity, % stable)
   - Error bars: SEM across neurons

2. **Line plot: Yield over sessions**
   - X: Session date
   - Y: % included
   - Separate lines by brain area
   - Trend line: Linear regression

---

## Part 3: Implementation Plan

### New Files to Create

1. **`code/plot_screening_metrics_session.m`**
   - Function: `plot_screening_metrics_session(session_data, screening_info, output_dir, session_name)`
   - Generates Figures 1.1-1.4 for a single session
   - Saves as `{session_name}_screening_metrics.pdf`

2. **`code/plot_screening_metrics_aggregate.m`**
   - Function: `plot_screening_metrics_aggregate(aggregate_data, output_dir)`
   - Generates Figures 2.1-2.6 for all sessions combined
   - Requires aggregated data structure (collect from all sessions)
   - Saves as `aggregate_screening_metrics.pdf`

3. **`code/aggregate_screening_data.m`**
   - Function: `aggregate_data = aggregate_screening_data(manifest, project_root)`
   - Loads screening_info from all sessions
   - Combines into single structure for aggregate plotting
   - Returns struct with fields for each metric across all neurons

### Modifications to Existing Files

1. **`code/run_4factors_analysis.m`**
   - Add toggle: `generate_screening_metrics_plots = true/false`
   - After screening step, call `plot_screening_metrics_session()` if enabled
   - After all sessions processed, call aggregate plotting functions

### Data Structure for Aggregation

```matlab
aggregate_data = struct(...
    'session_ids', {{}}, ...        % Cell array of session IDs
    'brain_areas', {{}}, ...        % Cell array ('SC' or 'SNc')
    'subjects', {{}}, ...           % Cell array of subject names
    'session_dates', [], ...        % Datetime array
    'cluster_ids', [], ...          % Global neuron identifiers
    'mean_firing_rates', [], ...    % [nNeurons x 1]
    'sparsity', [], ...             % [nNeurons x 1]
    'stability_scores', [], ...     % [nNeurons x 1]
    'presence_ratios', [], ...      % [nNeurons x 1]
    'stable_period_frac', [], ...   % [nNeurons x 1] fraction of trials
    'is_modulated', [], ...         % [nNeurons x 1] logical
    'is_included', [], ...          % [nNeurons x 1] logical
    'exclusion_reasons', {{}}, ...  % Cell array of strings
    'neuron_classes', {{}}, ...     % Cell array ('classic', 'interneuron', 'modulated', 'excluded')
    'fr_thresholds', [], ...        % [nNeurons x 1] threshold used for each neuron
    'session_indices', [] ...       % [nNeurons x 1] index into session_ids
);
```

### Integration with Pipeline

```matlab
% In run_4factors_analysis.m, after screening:
if generate_screening_metrics_plots
    plot_screening_metrics_session(session_data, screening_info, ...
        fullfile(project_root, 'figures'), session_id);
end

% After all sessions loop:
if generate_screening_metrics_plots
    aggregate_data = aggregate_screening_data(manifest, project_root);
    plot_screening_metrics_aggregate(aggregate_data, ...
        fullfile(project_root, 'figures'));
end
```

### Output Files

Session-level:
- `figures/{session_id}_screening_metrics.pdf` - 4-page PDF per session

Aggregate:
- `figures/aggregate_screening_metrics.pdf` - Multi-page PDF with all aggregate figures

---

## Implementation Order

1. **Phase 1:** Create `plot_screening_metrics_session.m` (Figures 1.1-1.4)
2. **Phase 2:** Create `aggregate_screening_data.m` (data collection utility)
3. **Phase 3:** Create `plot_screening_metrics_aggregate.m` (Figures 2.1-2.6)
4. **Phase 4:** Integrate into `run_4factors_analysis.m` with toggle
5. **Phase 5:** Test on all sessions and verify outputs

---

## Notes

- All figures use `exportgraphics` with `'Append', true` for multi-page PDFs
- Color scheme: Green = included, Red/Gray = excluded, Blue = threshold lines
- Log scale for FR axis when range spans multiple orders of magnitude
- Figure visibility set to 'off' for batch processing
- Consistent figure sizing: 1200x800 pixels for aggregate, 900x700 for session-level
