# **Project Progress and To-Do List**

This document tracks the development of the analysis pipeline. It will be updated as tasks are completed.

### **Part 0: Project Setup and Repository Initialization**

* \[x] **(USER)** Create a GitHub repository with a standard MATLAB .gitignore file.
* \[x] **(USER)** Create the initial directory structure: /code, /config, /docs, /output.
* \[x] **(USER)** Add existing helper functions (arrayROC.m, performCCA.m, barStairsFill.m, unitQualityAssesment\_script.m, makePSTHgrid.m) to the /code/utils/ directory.
* \[x] **(USER)** Place the initial README.md and PROGRESS\_AND\_TODO.md files in the repository root.
* \[x] **(USER)** Revise config/session\_manifest.csv to include new columns:
  * has\_gSac\_4factors (boolean/0-1)  
  * has\_gSac\_jph (boolean/0-1)  
  * neuron\_summary\_pdf\_status (e.g., 'pending', 'complete')  
  * behavioral\_analysis\_status (e.g., 'pending', 'complete')  
  * decoding\_analysis\_status (e.g., 'pending', 'complete')

### **Part 1: Data Loading & Preparation Module**

* \[x] **Task 1.1:** Create \+data\_handling/load\_session.m. This function will take a unique\_id string, read config/session\_manifest.csv, load the corresponding session\_data.mat file, and return both the session\_data struct and a struct of relevant metadata from the manifest row.
* \[x] **Task 1.2:** Create \+analysis/select\_neurons.m. This function acts as a controller to select neurons for analysis by calling area-specific screening functions. This is a more advanced implementation than the originally planned `get_good_neurons.m`.
* \[x] **Task 1.2a (New):** Create \+analysis/screen\_da\_neurons.m. Selects putative DA neurons based on firing rate and waveform characteristics.
* \[x] **Task 1.2b (New):** Create \+analysis/screen\_sc\_neurons.m. Implements an inclusive, multi-group method to identify task-modulated SC neurons.
* \[x] **Task 1.3:** Create \+analysis/get\_spike\_times.m. This function will take session\_data.spikes and a vector of cluster\_ids, and return a cell array where each cell contains the spike timestamps for one neuron. *(Note: Moved to +analysis package)*.
* \[x] **Task 1.4:** Create \+utils/alignAndBinSpikes.m. This function takes spike times, event alignment times, and a time window and returns a binned spike count matrix. *(Note: Implemented in utils package)*.

### **Part 2: Data Exploration & Quality Control**

* \[x] **Task 2.1:** The functionality for generating summary plots is partially implemented within the screening functions (`screen_da_neurons.m` and `screen_sc_neurons.m`), which create diagnostic plots for their respective analyses. A separate `generate_neuron_summary_pdfs.m` script does not exist.

### **Part 3: Behavioral Analysis Module**

* \[ \] **Task 3.1:** Create a script/function \+analysis/analyze\_behavior.m. This script will load data for a session, extract relevant variables from trialInfo (SRT, trialEndState, peakVel, stimType, reward, etc.), and use fitlme to assess the main effects and interactions of the task factors on behavior.

### **Part 4: Single-Neuron Analysis Module**

* \[ \] **Task 4.1:** Create \+analysis/run\_time\_resolved\_roc.m. This function will take binned spike data for two conditions and use utils/arrayROC.m to compute the time-resolved AUC.  
* \[ \] **Task 4.2:** Create \+analysis/run\_windowed\_roc.m. This function will average binned spike data over a specified time window and then compute a single AUC value per neuron using utils/arrayROC.m.

### **Part 5: Population Decoding Module**

* \[ \] **Task 5.1:** Create \+analysis/run\_time\_resolved\_svm.m. This function will take population data (nNeurons x nTrials x nTimeBins), and for each time bin, will train and evaluate a linear SVM classifier with k-fold cross-validation using fitcsvm.  
* \[ \] **Task 5.2:** Create \+analysis/run\_windowed\_svm.m. This will average population data over a time window and perform a single "snapshot" classification.

### **Part 6: Dimensionality Reduction Module**

* \[ \] **Task 6.1:** Create \+analysis/run\_dPCA.m. This function will format population data and apply the dPCA toolbox to identify and visualize components related to task factors (salience, reward, identity, etc.).

### **Part 7: Connectivity Analysis Module**

* \[ \] **Task 7.1:** Create \+analysis/run\_CCA.m. This function will take two population data matrices (e.g., SC and SNc activity), perform canonical correlation analysis using canoncorr, and return the results.