# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MATLAB-based neural data analysis pipeline investigating how Superior Colliculus (SC) and Substantia Nigra pars compacta (SNc) integrate cognitive/sensory signals (physical salience, reward value, spatial probability, target identity) into a "priority map" guiding saccadic eye movements. Processes standardized `session_data.mat` files from the upstream neuro-preprocessing-pipeline.

## Running the Pipeline

Execute MATLAB scripts via `-batch`:
```bash
matlab -batch "run('code/run_4factors_analysis.m')"       # Stage 1: Per-session analysis
matlab -batch "run('code/aggregate_analysis_results.m')"  # Stage 2: Data aggregation
matlab -batch "run('code/run_plotting_pipeline.m')"       # Stage 3: Generate figures
```

## Three-Stage Pipeline Architecture

```
session_data.mat → [Stage 1: Per-Session] → data/processed/{session_id}/
                                                      ↓
                   [Stage 2: Aggregation] → aggregated_analysis_data.mat
                                                      ↓
                   [Stage 3: Plotting]    → figures/
```

**Stage 1** (`run_4factors_analysis.m`): Loads sessions from `config/session_manifest.csv`, screens neurons, bins spikes, runs ROC/ANOVA/decoding analyses, saves per-session results.

**Stage 2** (`aggregate_analysis_results.m`): Concatenates all session results into plot-ready structures.

**Stage 3** (`run_plotting_pipeline.m`): Generates population summary figures.

## Single Source of Truth

`code/utils/define_task_conditions.m` defines the entire analysis plan:
- Alignment events: fixOn, targetOn, fixOff, saccadeOnset, reward
- Neuron inclusion thresholds (firing rate, sparsity, modulation)
- Two-model approach for confounded factors (Image Trials vs Bullseye Trials)
- Condition masks for all task variables

**Consult this file before modifying any analysis logic.**

## Key Conventions

### Path Setup (every script)
```matlab
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
```

### Coding Style
- Max 75 characters per line (use `...` for continuation)
- Comments above code, not to the right
- Use `%%` section markers for organization

### Coordinate System
- `trialInfo.targetTheta`: Stored as degrees × 10 — **divide by 10 before use**
- 0° = positive x-axis (3 o'clock), counter-clockwise increasing

### Neuron Classification

**SC neurons** (`screen_sc_neurons.m`):
- `classic`: Task-modulated + spatial selectivity
- `interneuron`: Visually suppressed or broadly activated
- `excluded`: Low FR, sparse, or not task-modulated

**SNc neurons** (`screen_snc_neurons.m`):
- `modulated`: Significant activity change + FR > 1 sp/s + not sparse
- `excluded`: Fails criteria

### SNc Subregion Labeling
Grid Y-coordinate from manifest:
- Y ∈ {3, 4} → `rvmSNc` (rostral-ventral-medial)
- Y ∈ {1, 2} → `cdlSNc` (caudal-dorsal-lateral)

### Plotting
- Use `barStairsFill.m` for PSTHs (not `plot`)
- Use `mySubPlot([a,b,c])` instead of `subplot(a,b,c)`
- Use `richColors.m` for consistent color palette
- Set `'Interpreter', 'none'` for titles with underscores

## Key Documentation

| Document | Purpose |
|----------|---------|
| `docs/design/analysis_rationale.md` | Scientific hypotheses and analysis logic |
| `docs/design/implementation_spec.md` | Data structures and function signatures |
| `docs/preprocessing_docs/session_data_dictionary.md` | Input `session_data.mat` structure |
| `docs/analysis_data_dictionary.md` | Analysis output structures |
| `AGENTS.md` | Extended AI assistant instructions |
| `PROGRESS_AND_TODO.md` | Development roadmap |

## Dependencies

- MATLAB R2021b+
- Statistics and Machine Learning Toolbox (required)
- Parallel Computing Toolbox (optional, for performance)
