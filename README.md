# **SC-SNc Saccade Priority Map Analysis Pipeline**

## **1\. Project Goal**

This repository contains a MATLAB-based analysis pipeline designed to investigate the neural mechanisms underlying saccadic eye movements. The primary scientific goal is to understand how the Superior Colliculus (SC) and Substantia Nigra pars compacta (SNc) interact to integrate various cognitive and sensory signals—such as physical salience, reward value, spatial probability, and target identity—into a unified "priority map" that guides behavior.

This pipeline processes standardized session\_data.mat files, which are the output of the neuro-preprocessing-pipeline.

## **2\. Repository Structure**

The project is organized into the following directories:

*   **/code/**: Contains all MATLAB source code for the analysis.
    *   **/utils/**: Contains general-purpose helper functions (e.g., `arrayROC.m`, `barStairsFill.m`).
    *   **/examples/**: Contains scripts that demonstrate how to use individual functions.
    *   All other analysis scripts are located in the root of the `/code/` directory.
*   **/config/**: Contains configuration files.
    *   `session_manifest.csv`: The central control file for the pipeline. It lists all available recording sessions and their metadata, and tracks the analysis progress.
*   **/docs/**: Contains project documentation, including data dictionaries.
*   **/figures/**: The destination for generated figures, such as per-neuron diagnostic plots. This directory is created automatically.
*   **/output/**: The destination for all generated files (figures, analysis results, reports). This directory is ignored by Git.

## **3. Core Dependencies**

*   **MATLAB (R2021b or newer recommended)**
*   **Required Toolboxes:**
    *   Statistics and Machine Learning Toolbox

## **4. Workflow Overview**

The main entry point for the analysis is the `run_4factors_analysis.m` script located in the `/code/` directory. This script automates the entire pipeline, from data loading to analysis and saving results.

The workflow is as follows:

1.  **Configuration**: The script is driven by the `config/session_manifest.csv` file. This file lists all sessions to be processed and tracks the status of each analysis stage (`screening`, `dataprep`, `analysis`).
2.  **Execution**: Run the `run_4factors_analysis.m` script from the MATLAB command line. The script will:
    a.  Load the session manifest.
    b.  Load an `analysis_plan` from `define_task_conditions.m`.
    c.  Iterate through each session in the manifest.
    d.  For each session, it performs the following steps, skipping any that are already marked as 'complete' in the manifest:
        i.  **Load Data**: Loads the `_session_data.mat` file.
        ii. **Neuron Screening**: Selects task-modulated neurons by calling `screen_sc_neurons.m` or `screen_da_neurons.m`.
        iii. **Diagnostic PDF Generation**: Creates per-neuron summary plots using `generate_neuron_summary_pdf.m`.
        iv. **Core Data Preparation**: Prepares the data for analysis by calling `prepare_core_data.m`.
        v.  **Run Analyses**: Executes a series of analyses as defined in the `analysis_plan`, including:
            *   `analyze_baseline_comparison.m`
            *   `analyze_roc_comparison.m`
            *   `analyze_anova.m`
    e.  Saves the updated `session_data` structure and the modified manifest.

You can use the `force_rerun` structure at the top of `run_4factors_analysis.m` to force specific stages of the pipeline to re-run, even if they are marked as complete.

## **5\. Data and Path Conventions**

This project uses a standardized two-directory system to separate code from data.

*   **Project Directory**: This is the repository you cloned from Git. It contains all the MATLAB code (`.m` files), configuration files, and documentation.
*   **Data Directory**: This is a separate directory located on OneDrive where the large `_session_data.mat` files are stored. It is kept separate from the project repository to avoid versioning large data files.

All scripts in this project are designed to be run from the root of the Project Directory. The code includes a helper function (`utils.findOneDrive`) to automatically locate the Data Directory on your system, ensuring that the analysis pipeline can find the necessary data files without manual path configuration.

## **6. Key Function Reference**

This section provides a brief overview of the key functions in the analysis pipeline. All functions are located in the `/code/` directory unless otherwise specified.

### **Main Script**

*   **`run_4factors_analysis.m`**: The main script that orchestrates the entire analysis pipeline.

### **Core Analysis Functions**

*   **`define_task_conditions.m`** (`/code/utils/`): Defines the `analysis_plan` that governs which analyses are run. It also generates logical masks to identify trials belonging to specific experimental conditions.
*   **`screen_sc_neurons.m`**: Identifies task-modulated Superior Colliculus (SC) neurons.
*   **`screen_da_neurons.m`**: Selects putative dopamine (DA) neurons based on electrophysiological criteria.
*   **`prepare_core_data.m`**: Aligns and bins spike data to create the `core_data` structure used by all subsequent analyses.
*   **`analyze_baseline_comparison.m`**: Compares post-event firing rates to a pre-event baseline.
*   **`analyze_roc_comparison.m`**: Performs ROC analysis to compare firing rates between two conditions over time.
*   **`analyze_anova.m`**: Performs N-way ANOVA on the data.
*   **`determine_rf_location.m`**: Determines the receptive field (RF) location for neurons.
*   **`generate_neuron_summary_pdf.m`**: Generates a multi-page PDF with diagnostic plots for each selected neuron.

### **Utility Functions (`/code/utils/`)**

*   **`alignAndBinSpikes.m`**: A core utility for creating binned spike count matrices.
*   **`arrayROC.m`**: Calculates the Receiver Operating Characteristic (ROC) and Area Under the Curve (AUC) for time-resolved data.
*   **`calculate_baseline_fr.m`**: Computes the baseline firing rate for each neuron.
*   **`findOneDrive.m`**: A helper function that automatically finds the path to the user's OneDrive directory.
*   **`initCodes.m`**: Initializes a structure with task-related codes.