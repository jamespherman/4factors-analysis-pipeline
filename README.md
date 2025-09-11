# **SC-SNc Saccade Priority Map Analysis Pipeline**

## **1\. Project Goal**

This repository contains a MATLAB-based analysis pipeline designed to investigate the neural mechanisms underlying saccadic eye movements. The primary scientific goal is to understand how the Superior Colliculus (SC) and Substantia Nigra pars compacta (SNc) interact to integrate various cognitive and sensory signals—such as physical salience, reward value, spatial probability, and target identity—into a unified "priority map" that guides behavior.

This pipeline processes standardized session\_data.mat files, which are the output of the neuro-preprocessing-pipeline.

## **2\. Repository Structure**

The project is organized into the following directories:

* **/code/**: Contains all MATLAB source code for the analysis.  
  * main\_analysis\_pipeline.m: The master script to run the entire analysis pipeline.  
  * \+data\_handling/: Package for functions related to loading, preparing, and organizing data.  
  * \+analysis/: Package containing modules for each major analysis type (behavioral, single-neuron ROC, population decoding, CCA, etc.).  
  * \+visualization/: Package for all plotting and figure generation functions.  
  * utils/: Directory for general-purpose helper functions (e.g., arrayROC.m, barStairsFill.m).  
* **/config/**: Contains configuration files.  
  * session\_manifest.csv: The central control file for the pipeline. It lists all available recording sessions and their metadata, and is used to select data for analysis.  
* **/docs/**: Contains project documentation, including the master development plan.  
* **/output/**: The destination for all generated files (figures, analysis results, reports). This directory is ignored by Git.

## **3\. Core Dependencies**

* **MATLAB (R2021b or newer recommended)**  
* **Required Toolboxes:**  
  * Statistics and Machine Learning Toolbox  
* **External Toolboxes:**  
  * [dPCA Toolbox](https://github.com/machenslab/dPCA) (to be added)

## **4\. Workflow Overview**

The analysis pipeline is controlled by the session\_manifest.csv file and executed via the main script.

1. **Configure Manifest**: Before running, modify config/session\_manifest.csv to select the desired sessions for analysis. The manifest should be updated to include boolean flags for the presence of specific tasks (e.g., has\_gSac\_4factors) and status columns for tracking analysis progress (e.g., behavioral\_analysis\_status).  
2. **Run Pipeline**: Execute the master script from the MATLAB command window:  
   \>\> main\_analysis\_pipeline

3. **Review Output**: The script will iterate through the selected sessions, perform the specified analyses, and save all outputs (e.g., per-neuron PDF summaries, statistical results, figures) to the /output directory, organized by session unique\_id.

## **5\. Data and Path Conventions**

This project uses a standardized two-directory system to separate code from data.

*   **Project Directory**: This is the repository you cloned from Git. It contains all the MATLAB code (`.m` files), configuration files, and documentation.
*   **Data Directory**: This is a separate directory located on OneDrive where the large `_session_data.mat` files are stored. It is kept separate from the project repository to avoid versioning large data files.

All scripts in this project are designed to be run from the root of the Project Directory. The code includes a helper function (`utils.findOneDrive`) to automatically locate the Data Directory on your system, ensuring that the analysis pipeline can find the necessary data files without manual path configuration.