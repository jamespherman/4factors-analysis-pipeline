# **Project Progress and To-Do List**

This document tracks the development of the analysis pipeline. It has been updated to reflect the current state of the codebase, which is primarily driven by the `run_4factors_analysis.m` script.

### **Part 0: Project Setup and Repository Initialization**

*   [x] **(USER)** Create a GitHub repository with a standard MATLAB .gitignore file.
*   [x] **(USER)** Create the initial directory structure: /code, /config, /docs, /output.
*   [x] **(USER)** Add existing helper functions to the `/code/utils/` directory.
*   [x] **(USER)** Place the initial README.md and PROGRESS_AND_TODO.md files in the repository root.
*   [x] **(USER)** Revise `config/session_manifest.csv` to include status columns (`screening_status`, `dataprep_status`, `analysis_status`).

### **Part 1: Core Pipeline Implementation**

*   [x] **Task 1.1: Master Script (`run_4factors_analysis.m`)**
    *   Reads `config/session_manifest.csv`.
    *   Loads a dynamic `analysis_plan` from `define_task_conditions.m`.
    *   Iterates through sessions, checking the status of each step.
    *   Saves updated data and manifest.
*   [x] **Task 1.2: Neuron Screening (`screen_da_neurons.m`, `screen_sc_neurons.m`)**
    *   Selects putative DA neurons based on firing rate and waveform.
    *   Identifies task-modulated SC neurons using a multi-group method.
*   [x] **Task 1.3: Data Preparation (`prepare_core_data.m`)**
    *   Aligns and bins spike data for all required events and conditions, creating the `core_data` structure.
*   [x] **Task 1.4: Analysis Plan (`define_task_conditions.m`)**
    *   When called with no arguments, returns the `analysis_plan` struct.
    *   When called with `session_data`, returns a struct of logical masks for all experimental conditions.

### **Part 2: Single-Neuron Analysis and Diagnostics**

*   [x] **Task 2.1: Per-Neuron Diagnostic PDFs (`generate_neuron_summary_pdf.m`)**
    *   Generates a multi-page PDF with key diagnostic plots (PSTHs, etc.) for each selected neuron in a session.
*   [x] **Task 2.2: Baseline vs. Activity Comparison (`analyze_baseline_comparison.m`)**
    *   Performs ROC analysis to compare post-event firing rates to a pre-event baseline period for various conditions.
*   [x] **Task 2.3: Time-Resolved ROC Analysis (`analyze_roc_comparison.m`)**
    *   Compares firing rates between two experimental conditions over time using a sliding window ROC analysis (`utils/arrayROC.m`).
*   [x] **Task 2.4: N-way ANOVA (`analyze_anova.m`)**
    *   Performs an N-way ANOVA to investigate main effects and interactions of task factors on neural activity.

### **Part 3: Aggregation and Plotting**

*   [x] **Task 3.1: Aggregate Analysis Results (`aggregate_analysis_results.m`)**
    *   Iterates through all completed session data files.
    *   Extracts and concatenates the results from all single-session analyses into a single, aggregated data file (`aggregated_analysis_results.mat`).
*   [x] **Task 3.2: Plot Aggregated Results**
    *   `plot_aggregated_baseline_comparison.m`
    *   `plot_aggregated_roc_comparison.m`
    *   `plot_aggregated_anova.m`
    *   These scripts load the aggregated data and generate summary plots across all sessions.

### **Part 4: Future Work and Optimization (To-Do)**

*   [ ] **Task 4.1: Behavioral Analysis Module (`analyze_behavior.m`)**
    *   Create a script to analyze behavioral data (e.g., reaction times, success rates) and assess the effects of task factors using `fitlme`.
*   [ ] **Task 4.2: Population Decoding Module**
    *   Implement functions to perform time-resolved and windowed SVM classification to decode task variables from neural population activity.
*   [ ] **Task 4.3: Dimensionality Reduction Module**
    *   Implement dPCA or other dimensionality reduction techniques to identify and visualize population-level coding dimensions related to task factors.
*   [ ] **Task 4.4: Code Refactoring and Modularization**
    *   The current structure has many scripts in the root `/code/` directory. Consider refactoring into a more modular package structure (e.g., `+analysis`, `+plotting`) to improve organization and maintainability, as was originally envisioned.
*   [ ] **Task 4.5: Update Function-Level Documentation**
    *   Ensure all functions have complete and accurate header comments explaining their purpose, inputs, outputs, and any important assumptions.

### **run_4factors_analysis.m Audit**

*This section documents findings from a code audit performed on `run_4factors_analysis.m`. The goal was to identify bugs, inconsistencies, and logical errors without implementing fixes.*

1.  **Inconsistent Project Root Path:**
    *   **File:** `code/run_4factors_analysis.m`
    *   **Line:** 28
    *   **Issue:** The script defines `project_root` using `fileparts`, which contradicts the `findOneDrive`-based method mandated in `code/AGENTS.md`. This could lead to pathing errors if the execution environment changes.

2.  **Improper MATLAB Path Addition:**
    *   **File:** `code/run_4factors_analysis.m`
    *   **Line:** 29
    *   **Issue:** The script adds the entire project root directory to the MATLAB path. This is poor practice as it can lead to function shadowing and conflicts with built-in MATLAB functions.

3.  **Access to Non-Existent Field:**
    *   **File:** `code/run_4factors_analysis.m`
    *   **Line:** 88
    *   **Issue:** The script checks for `session_data.metrics`, a field that is not defined in the `docs/preprocessing_docs/session_data_dictionary.md`. This may be a remnant of an older data structure and will always be false.

4.  **Hardcoded Event Name in Idempotency Check:**
    *   **File:** `code/run_4factors_analysis.m`
    *   **Lines:** 400, 513
    *   **Issue:** The idempotency check for the "Baseline Comparison" analysis and the final manifest verification hardcode the event name `'targetOn'`. This makes the logic brittle and will cause unnecessary re-computation or failed completion marking if a session's analysis plan does not involve that specific event.

5.  **Flawed "Dry Run" and Execution Logic:**
    *   **File:** `code/run_4factors_analysis.m`
    *   **Lines:** 145-285
    *   **Issue:** The population decoding execution block (lines 145-258) is incorrectly placed within the "dry run" section that calculates the total number of steps. This means progress reporting (`Step X of Y`) will be inaccurate, as the total `Y` is not finalized before the steps begin.

6.  **Uninitialized Variable Error:**
    *   **File:** `code/run_4factors_analysis.m`
    *   **Lines:** 203, 238, 286
    *   **Issue:** The variable `step_counter` is used inside the misplaced population decoding block (e.g., lines 203, 238) *before* it is initialized on line 286. This will cause a fatal script error.

7.  **Inconsistent Idempotency Check Implementation:**
    *   **File:** `code/run_4factors_analysis.m`
    *   **Lines:** 396-484
    *   **Issue:** The logic for skipping completed analyses is inconsistent. "Baseline Comparison" uses a hardcoded event, "ROC Comparison" correctly uses the event from the analysis plan, and "ANOVA" uses a simple `isfield` check. This inconsistency makes the script harder to maintain and debug.

8.  **Flawed Population Decoding Idempotency:**
    *   **File:** `code/run_4factors_analysis.m`
    *   **Lines:** 195-224
    *   **Issue:** The check to skip model training is insufficient. It reruns training if a previous attempt completed training but failed during the subsequent testing phase, which is inefficient. The logic does not properly distinguish between a run that has not started and one that is partially complete.