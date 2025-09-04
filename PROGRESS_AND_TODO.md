# **Project Progress and To-Do List**

This document tracks the development of the analysis pipeline. It will be updated as tasks are completed.

### **Part 0: Project Setup and Repository Initialization**

* \[ \] **(USER)** Create a GitHub repository with a standard MATLAB .gitignore file.  
* \[ \] **(USER)** Create the initial directory structure: /code, /config, /docs, /output.  
* \[ \] **(USER)** Add existing helper functions (arrayROC.m, performCCA.m, barStairsFill.m, unitQualityAssesment\_script.m, makePSTHgrid.m) to the /code/utils/ directory.  
* \[ \] **(USER)** Place the initial README.md and PROGRESS\_AND\_TODO.md files in the repository root.  
* \[ \] **(USER)** Revise config/session\_manifest.csv to include new columns:  
  * has\_gSac\_4factors (boolean/0-1)  
  * has\_gSac\_jph (boolean/0-1)  
  * neuron\_summary\_pdf\_status (e.g., 'pending', 'complete')  
  * behavioral\_analysis\_status (e.g., 'pending', 'complete')  
  * decoding\_analysis\_status (e.g., 'pending', 'complete')

### **Part 1: Data Loading & Preparation Module**

* \[ \] **Task 1.1:** Create \+data\_handling/load\_session.m. This function will take a unique\_id string, read config/session\_manifest.csv, load the corresponding session\_data.mat file, and return both the session\_data struct and a struct of relevant metadata from the manifest row.  
* \[ \] **Task 1.2:** Create \+utils/get\_good\_neurons.m. This function will take a session\_data.spikes struct and return a vector of cluster\_ids for neurons marked as 'good'.  
* \[ \] **Task 1.3:** Create \+data\_handling/get\_spike\_times.m. This function will take session\_data.spikes and a vector of cluster\_ids, and return a cell array where each cell contains the spike timestamps for one neuron.  
* \[ \] **Task 1.4:** Create \+data\_handling/bin\_spike\_data.m. This function will take spike times, event alignment times, and a time window (e.g., \[-0.5, 1.0\]), and return a binned spike count matrix (nTrials x nTimeBins) for a single neuron.

### **Part 2: Data Exploration & Quality Control**

* \[ \] **Task 2.1:** Create a script/function generate\_neuron\_summary\_pdfs.m. For each 'good' neuron in a session, this script will generate and save a single-page PDF containing:  
  * Mean waveform on the peak channel.  
  * Inter-Spike Interval (ISI) histogram.  
  * A grid of PSTHs for key trial conditions from the gSac\_4factors task, using subplot (not TiledLayout) and barStairsFill.m for plotting.

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