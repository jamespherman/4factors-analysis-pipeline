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

The analysis pipeline is executed as a sequence of three main scripts. This multi-step design separates per-session processing from data aggregation and plotting, providing a clear and modular workflow.

The basic workflow is as follows:

1.  **`run_4factors_analysis.m`**: This is the first stage of the pipeline. It iterates through each session listed in `config/session_manifest.csv` and performs all the core analyses on a per-session basis. For each session, it:
    *   Loads the raw `_session_data.mat` file.
    *   Screens neurons to identify task-modulated units.
    *   Prepares the data for analysis by creating a `core_data` structure.
    *   Executes a series of analyses as defined in `define_task_conditions.m`, including ROC, ANOVA, and baseline comparisons.
    *   Saves the results back into a processed `_session_data.mat` file in the `/output/` directory and updates the session manifest.

2.  **`aggregate_analysis_results.m`**: This script runs after all sessions have been processed by `run_4factors_analysis.m`. It loops through the processed `_session_data.mat` files and compiles the results from all sessions into a single, aggregated data file: `aggregated_analysis_data.mat`. This file is structured to be "plot-ready" for the final stage.

3.  **`run_plotting_pipeline.m`**: This is the final step. It loads the `aggregated_analysis_data.mat` file and uses the aggregated data to generate all the summary figures for the project.

This sequential process ensures that data is first processed on a granular, per-session level, then efficiently aggregated, and finally plotted.

## **5\. Data and Path Conventions**

This project uses a standardized two-directory system to separate code from data.

*   **Project Directory**: This is the repository you cloned from Git. It contains all the MATLAB code (`.m` files), configuration files, and documentation.
*   **Data Directory**: This is a separate directory located on OneDrive where the large `_session_data.mat` files are stored. It is kept separate from the project repository to avoid versioning large data files.

All scripts in this project are designed to be run from the root of the Project Directory. The code includes a helper function (`utils.findOneDrive`) to automatically locate the Data Directory on your system, ensuring that the analysis pipeline can find the necessary data files without manual path configuration.

## **6. Key Function Reference**

This section provides a brief overview of the key functions in the analysis pipeline. All functions are located in the `/code/` directory unless otherwise specified.

### **Main Pipeline Scripts**

*   **`run_4factors_analysis.m`**: The first stage of the pipeline. It runs all analyses on a per-session basis.
*   **`aggregate_analysis_results.m`**: The second stage. It aggregates the results from all individual sessions into a single data file.
*   **`run_plotting_pipeline.m`**: The final stage. It generates all summary figures from the aggregated data.

### **Core Analysis Functions**

*   **`define_task_conditions.m`** (`/code/utils/`): Defines the `analysis_plan` that governs which analyses are run. It also generates logical masks to identify trials belonging to specific experimental conditions.
*   **`screen_sc_neurons.m`**: Identifies task-modulated Superior Colliculus (SC) neurons.
*   **`screen_da_neurons.m`**: Selects putative dopamine (DA) neurons based on electrophysiological criteria.
*   **`prepare_core_data.m`**: Aligns and bins spike data to create the `core_data` structure used by all subsequent analyses.
*   **`analyze_baseline_comparison.m`**: Compares post-event firing rates to a pre-event baseline.
*   **`analyze_roc_comparison.m`**: Performs ROC analysis to compare firing rates between two conditions over time.
*   **`analyze_anova.m`**: Performs N-way ANOVA on the data.
*   **`generate_neuron_summary_pdf.m`**: Generates a multi-page PDF with diagnostic plots for each selected neuron.

### **Utility Functions (`/code/utils/`)**

*   **`alignAndBinSpikes.m`**: A core utility for creating binned spike count matrices.
*   **`arrayROC.m`**: Calculates the Receiver Operating Characteristic (ROC) and Area Under the Curve (AUC) for time-resolved data.
*   **`calculate_baseline_fr.m`**: Computes the mean firing rate for each neuron over a specified window for a given set of trials.
*   **`findOneDrive.m`**: A helper function that automatically finds the path to the user's OneDrive directory.
*   **`initCodes.m`**: Initializes a structure with task-related codes.