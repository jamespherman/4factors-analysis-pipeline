# **SC-SNc Saccade Priority Map Analysis Pipeline**

## **1\. Project Goal**

This repository contains a MATLAB-based analysis pipeline designed to investigate the neural mechanisms underlying saccadic eye movements. The primary scientific goal is to understand how the Superior Colliculus (SC) and Substantia Nigra pars compacta (SNc) interact to integrate various cognitive and sensory signals—such as physical salience, reward value, spatial probability, and target identity—into a unified "priority map" that guides behavior.

This pipeline processes standardized session\_data.mat files, which are the output of the neuro-preprocessing-pipeline.

## **2\. Repository Structure**

The project is organized into the following directories:

* **/code/**: Contains all MATLAB source code for the analysis.  
  * \+data\_handling/: Package for functions related to loading and preparing data.
  * \+analysis/: Package containing modules for core analysis tasks like neuron screening.
  * \+utils/: Directory for general-purpose helper functions (e.g., arrayROC.m, barStairsFill.m).
  * \+examples/: Contains scripts that demonstrate how to use the core functions to perform analyses.
* **/config/**: Contains configuration files.  
  * session\_manifest.csv: The central control file for the pipeline. It lists all available recording sessions and their metadata.
* **/docs/**: Contains project documentation.
* **/output/**: The destination for all generated files (figures, analysis results, reports). This directory is ignored by Git.

## **3\. Core Dependencies**

* **MATLAB (R2021b or newer recommended)**  
* **Required Toolboxes:**  
  * Statistics and Machine Learning Toolbox

## **4\. Workflow Overview**

This repository does not currently have a single master script to run all analyses. Instead, the analysis is performed by running scripts from the **/code/examples/** directory. These scripts demonstrate the intended workflow:

1. **Load Data**: Use `data_handling.load_session` to load a specific session's data by referencing its `unique_id` from `config/session_manifest.csv`.
2. **Select Neurons**: Use `analysis.select_neurons` to identify task-modulated neurons based on the recorded brain area (SC or SNc). This function calls area-specific screening functions (`screen_sc_neurons` or `screen_da_neurons`).
3. **Perform Analysis**: Use the various functions and example scripts to perform specific analyses, such as calculating firing rates, plotting PSTHs, or computing ROC curves.

Please see the scripts in the `/code/examples/` directory for concrete examples of how to use the functions in this pipeline.

## **5\. Data and Path Conventions**

This project uses a standardized two-directory system to separate code from data.

*   **Project Directory**: This is the repository you cloned from Git. It contains all the MATLAB code (`.m` files), configuration files, and documentation.
*   **Data Directory**: This is a separate directory located on OneDrive where the large `_session_data.mat` files are stored. It is kept separate from the project repository to avoid versioning large data files.

All scripts in this project are designed to be run from the root of the Project Directory. The code includes a helper function (`utils.findOneDrive`) to automatically locate the Data Directory on your system, ensuring that the analysis pipeline can find the necessary data files without manual path configuration.

## **6\. Function Reference**

This section provides a brief overview of the key functions in the analysis pipeline.

### **`+data_handling`**

*   **`load_session.m`**: Loads a session's `_session_data.mat` file and merges it with the corresponding metadata from `session_manifest.csv`.

### **`+analysis`**

*   **`select_neurons.m`**: A top-level function that selects neurons for analysis based on the brain area recorded. It calls one of the area-specific screening functions below.
*   **`screen_sc_neurons.m`**: Implements an inclusive, multi-group method to identify task-modulated Superior Colliculus (SC) neurons. It identifies neurons that show significant firing rate changes during different task epochs (e.g., visual, delay, saccade) across various trial conditions.
*   **`screen_da_neurons.m`**: Selects putative dopamine (DA) neurons in the Substantia Nigra pars compacta (SNc) based on electrophysiological criteria, specifically a low baseline firing rate (< 20 sp/s) and a broad waveform.
*   **`get_spike_times.m`**: Extracts and returns a cell array of spike timestamps for a given list of neuron IDs.
*   **`define_task_conditions.m`**: Defines logical masks to identify trials belonging to specific experimental conditions (e.g., based on target location, reward value).
*   **`determine_rf_location.m`**: Determines the receptive field (RF) location for neurons, typically by comparing activity for contralateral vs. ipsilateral targets.

### **`+utils` (Key Functions)**

*   **`alignAndBinSpikes.m`**: A core utility for creating spike count matrices. It takes spike times, alignment event times, and a time window, and returns a binned matrix of spike counts, which is fundamental for many analyses.
*   **`arrayROC.m`**: Calculates the Receiver Operating Characteristic (ROC) curve and the Area Under the Curve (AUC) for time-resolved data, allowing for analysis of how well a neuron's firing rate discriminates between two conditions over time.
*   **`calculate_baseline_fr.m`**: Computes the baseline firing rate for each neuron.
*   **`calculate_waveform_metrics.m`**: Calculates metrics from a neuron's average waveform shape, such as peak-to-trough duration.
*   **`findOneDrive.m`**: A helper function that automatically finds the path to the user's OneDrive directory, allowing for seamless access to the data directory.