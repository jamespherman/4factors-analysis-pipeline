# Agent Instructions

This document provides instructions for AI coding assistants (Claude Code, etc.) working on the 4factors analysis pipeline.

---

## Project Overview

This is a MATLAB-based neural data analysis pipeline for studying how stimulus-driven and goal-directed factors are encoded in primate Superior Colliculus (SC) and Substantia Nigra pars compacta (SNc). The pipeline analyzes data from the `gSac_4factors` memory-guided saccade task.

**For scientific context and analysis rationale, see:** `docs/design/analysis_rationale.md`

---

## Execution Environment

MATLAB code can be executed using the `-batch` option:
```bash
matlab -batch "run('script_name.m')"
```

Or for single expressions:
```bash
matlab -batch "disp(pwd); disp(version)"
```

Notes:
- `-batch` runs without the GUI and exits on completion
- Errors will return a non-zero exit code
- Working directory and path setup may need explicit handling
---

## Key Documentation

| Document | Location | Purpose |
|----------|----------|---------|
| Analysis Rationale | `docs/design/analysis_rationale.md` | Scientific hypotheses, analysis logic, interpretation |
| Implementation Spec | `docs/design/implementation_spec.md` | Data structures, function signatures, code patterns |
| Session Data Dictionary | `docs/preprocessing_docs/session_data_dictionary.md` | Structure of `session_data.mat` files |
| Analysis Data Dictionary | `docs/analysis_data_dictionary.md` | Structure of analysis outputs |
| Plotting Requirements | `docs/plotting_requirements.md` | Figure specifications |

**Consult these documents before modifying code that interacts with the relevant data structures.**

---

## Refactoring and Impact Analysis Protocol

The analysis pipeline follows a strict sequential order:

```
Define Conditions → Run Per-Session Analyses → Aggregate Results → Plot Results
```

**When modifying code, assess impact on all subsequent stages.** A change to `define_task_conditions.m` affects all stages; a change to a plotting function affects none.

Before writing any code, use this checklist:

1. **Definition**: Which files define the relevant data structures or configurations?
   - `define_task_conditions.m`, `.md` files in `/docs`

2. **Creation**: Which scripts create the data on a per-session basis?
   - `analyze_*.m` functions

3. **Orchestration**: Which scripts execute the creation process?
   - `run_4factors_analysis.m`

4. **Aggregation**: Which scripts pool data from multiple sessions?
   - `aggregate_analysis_results.m`

5. **Consumption**: Which scripts consume the final aggregated data?
   - `plot_aggregated_*.m` functions

---

## Data Structure Hierarchy of Truth

**Before writing or modifying any code that interacts with `session_data.mat` files, consult the data dictionaries.**

For analyses driven by `define_task_conditions.m`:

1. The `condition_defs` structure produced by `define_task_conditions.m` is the **absolute source of truth** for the analysis plan.
2. The `analysis_plan` definition in `docs/analysis_data_dictionary.md` serves as documentation.

**If you find discrepancies between code output and documentation, the code takes precedence.**

Key structure notes:
- Number of neurons: derived from `session_data.spikes.cluster_info`, not a `nClusters` field
- Event times: in `session_data.eventTimes` struct, not a `trials` struct
- Spike times: single vector (`spikes.times`) mapped via `spikes.clusters`, not cell array per neuron

---

## Coding Style

### Line Length
Maximum 75 characters per line. Use `...` for line continuation.

### Commenting
Place comments on the line **above** the code they describe, not to the right.

### Path Setup
Add the `utils` directory to the path at the beginning of scripts:

```matlab
%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
```

### Project Root
Access project root using:

```matlab
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
```

---

## File Structure

All new MATLAB scripts should begin with a standard header:

```matlab
%% my_script.m
%
% A brief one-line description of what the script does.
%
% A more detailed description can follow, explaining the inputs, outputs,
% and any key assumptions or dependencies.
%
% Author: Your Name
% Date: YYYY-MM-DD
%
```

Organize code into logical sections using `%%` headers:

```matlab
%% Setup
% ... setup code here ...

%% Main Analysis
% ... main analysis code here ...

%% Plotting
% ... plotting code here ...
```

---

## Analysis Conventions

### Selecting Trials by Task

Filter trials based on task using `session_data.trialInfo.taskCode`:

```matlab
codes = initCodes();
fourfactors_trials = session_data.trialInfo.taskCode == ...
    codes.uniqueTaskCode_4factors;
```

### Condition Mask Compatibility

Masks from `define_task_conditions.m` are sized for filtered trial sets. They can be directly applied to `core_data` arrays, which are also filtered.

### Two-Model Approach

The 4factors task has confounded factors requiring separate analysis models:
- **Image Trials Model**: Factors are reward, probability, identity
- **Bullseye Trials Model**: Factors are reward, probability, salience

See `docs/design/analysis_rationale.md` Section 1.3 for details.

### SNc Subregion Labeling

SNc neurons are labeled by subregion based on grid coordinates (from manifest):
- Y ∈ {3, 4} → `rvmSNc` (rostral-ventral-medial)
- Y ∈ {1, 2} → `cdlSNc` (caudal-dorsal-lateral)

**This labeling only applies when `brain_area == 'SNc'`.**

### gSac_jph Task Properties

Memory-guided saccade trials in `gSac_jph` are placed at the neuron's RF center:
1. `scSide` can be inferred from target location
2. All trials are effectively "contralateral"

---

## Neuron Screening

### SC Neurons (`screen_sc_neurons.m`)

```matlab
[selected_neurons, sig_epoch_comp, scSide] = screen_sc_neurons(session_data);
```

- Automatically determines recorded SC side
- Uses memory-guided saccade trials from `gSac_jph`
- Returns `sig_epoch_comp`: [nClusters × 3] boolean matrix for Visual/Delay/Saccade vs Baseline

### Neuron Classification Categories

| Category | Criteria |
|----------|----------|
| `visually_excited` | Significant increase in visual epoch (contralateral) |
| `visually_suppressed` | Significant decrease in visual epoch (contralateral) |
| `delay_modulated` | Significant change in delay epoch |
| `saccade_modulated` | Significant change in saccade epoch |
| `unmodulated` | No significant modulation |

---

## Plotting Conventions

### PSTHs
Use `barStairsFill.m` instead of `plot`:

```matlab
% Single PSTH
hBS = barStairsFill(time_vector, zeros(size(mean_psth)), mean_psth);
delete(hBS(2))
set(hBS(1), 'FaceColor', 'k');
set(hBS(3), 'Color', 'k')

% Multiple PSTHs on same axes
hold on;
hBS1 = barStairsFill(time_vector, zeros(size(psth1)), psth1);
delete(hBS1(1:2))
set(hBS1(3), 'Color', 'b');
```

### Multi-Panel Figures
- Only bottom row has X-axis tick labels
- Only left column has Y-axis tick labels
- Use `mySubPlot([a,b,c])` instead of `subplot(a,b,c)`
- Column widths proportional to time window duration

### Titles
- Use `sgtitle()` for main figure title
- Set `'Interpreter', 'none'` for titles with underscores/file paths

### Colors
Use `richColors.m` for consistent color palette (16 distinct colors).

---

## Coordinate System

- **`trialInfo.targetTheta`**: Stored as degrees × 10; **divide by 10 before use**
- **0°**: Positive x-axis (3 o'clock)
- **Direction**: Counter-clockwise increasing
- **Right Visual Field**: θ ∈ [0°, 90°) ∪ (270°, 360°)
- **Left Visual Field**: θ ∈ (90°, 270°)

---

## Current Implementation Priorities

See `PROGRESS_AND_TODO.md` for detailed task list. Current focus areas:

1. **SNc subregion labeling** - Add to condition masks
2. **Window-based ROC** - Fixed-window summary metrics
3. **Neuron metrics table** - For population scatter plots
4. **Scatter plots with tooltips** - Neuron identification
5. **Single-page neuron summary PDFs** - Comprehensive diagnostics
6. **Decoding plot reorganization** - Matrix layout for cross-factor tests
