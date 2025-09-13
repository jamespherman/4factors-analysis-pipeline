# Analysis Data Dictionary

This document defines the canonical data structure for the `session_data.analysis` field. All per-session analysis functions must save their results in this format.

The guiding principle is to organize data by **Analysis Type -> Alignment Event -> Condition/Comparison Name**.

## Top-Level Structure

-   `session_data.analysis`
    -   `.roc_comparison`: Results from `analyze_roc_comparison.m`
    -   `.baseline_comparison`: Results from `analyze_baseline_comparison.m`
    -   `.anova_results`: Results from `analyze_anova.m`
    -   `.selected_neurons`: List of neurons selected for analysis.
    -   `.core_data`: Processed data used for analyses.
    -   `.scSide` (for SC sessions): The determined side of the Superior Colliculus recording ('left' or 'right').
    -   `.sig_epoch_comparison` (for SC sessions): Significance of epoch comparisons from `screen_sc_neurons.m`.

---

## 1. ROC Comparison

Compares firing rates between two conditions.

-   **Path:** `session_data.analysis.roc_comparison.(eventName).(compName)`
-   **`eventName`**: `char` - The name of the alignment event (e.g., `'CUE_ON'`, `'outcomeOn'`).
-   **`compName`**: `char` - The name of the specific comparison (e.g., `'Dist_at_Cue'`, `'RPE_at_Outcome'`).

### Fields

-   `.sig`: `[nNeurons x nTimeBins] double` - Matrix of p-values from the ROC analysis at each time bin.
-   `.time_vector`: `[1 x nTimeBins] double` - The time vector corresponding to the columns of `.sig`.
-   `.cond_names`: `[1 x 2] cell` - Cell array containing the names of the two conditions that were compared (e.g., `{'is_normal_dist', 'is_uniform_dist'}`).

---

## 2. Baseline Comparison

Compares post-event firing rates to a pre-event baseline period.

-   **Path:** `session_data.analysis.baseline_comparison.(eventName).(condName)`
-   **`eventName`**: `char` - The name of the alignment event (e.g., `'CUE_ON'`, `'outcomeOn'`).
-   **`condName`**: `char` - The name of the condition being analyzed (e.g., `'is_normal_dist'`, `'is_common_reward_no_spe'`).

### Fields

-   `.sig`: `[nNeurons x nTimeBins] double` - Matrix of p-values from the ROC analysis comparing each post-event bin to the pooled baseline.
-   `.time_vector`: `[1 x nTimeBins] double` - The time vector for the post-event period, corresponding to the columns of `.sig`.

---

## 3. ANOVA Results

Results from the N-way ANOVA analysis.

-   **Path:** `session_data.analysis.anova_results.(eventName).(p_value_name)`
-   **`eventName`**: `char` - The name of the alignment event (e.g., `'CUE_ON'`).
-   **`p_value_name`**: `char` - The name of the factor or interaction term for which p-values are reported (e.g., `'p_interaction'`, `'p_value'`).

### Fields

The fields under `(p_value_name)` are typically `[nNeurons x nTimeBins]` double matrices containing the p-values for that factor over time.
