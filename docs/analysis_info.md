# Guide for the Data Analyst

This document provides essential information for analyzing the preprocessed data in this project. The final output of the pipeline is a `session_data.mat` file for each recording session, but this file does not contain all the necessary information for a complete analysis.

## Session Metadata: `config/session_manifest.csv`

The primary source of metadata for each session is the `config/session_manifest.csv` file. This file is the control center for the entire preprocessing pipeline. Each row in this file corresponds to a unique recording from a single probe and contains the following critical information:

*   **`unique_id`**: A unique identifier for the session (e.g., `Feynman_08_12_2025_SNc`). This ID is used to name the output files, including the `session_data.mat` file.
*   **`probe_type`**: The type of neural probe used.
*   **`brain_area`**: The targeted brain region.
*   **`channel_ordering`**: The physical layout of channels on the probe.
*   **`notes`**: A free-text field that often contains information about the specific behavioral task run during the session.

To begin an analysis, you should always start by consulting the manifest to identify the sessions of interest and retrieve their associated metadata.

## Locating Session Data Files

All processed data, including the final `_session_data.mat` file, is stored in a structured directory within your OneDrive. The root for all processed data is a folder named `Neuronal Data Analysis`. Inside this folder, there is a separate directory for each session, named with the session's `unique_id`.

The full path to a session's data folder is:
`{Your_OneDrive_Path}/Neuronal Data Analysis/{unique_id}/`

The pipeline includes a helper function, `utils.findOneDrive.m`, to get the root path of your OneDrive folder programmatically. You can combine this with a `unique_id` from the manifest to construct the full path to the data file you want to analyze.

Example (in MATLAB):
```matlab
% 1. Choose a unique_id from the session_manifest.csv
session_id = 'Feynman_08_12_2025_SNc';

% 2. Get the root of your OneDrive directory
oneDrivePath = utils.findOneDrive();

% 3. Construct the path to the session's data folder
sessionDataFolder = fullfile(oneDrivePath, 'Neuronal Data Analysis', session_id);

% 4. Construct the full path to the session_data.mat file
sessionDataFile = fullfile(sessionDataFolder, [session_id, '_session_data.mat']);

% 5. Check if the file exists and load it
if isfile(sessionDataFile)
    fprintf('Loading data from: %s\n', sessionDataFile);
    load(sessionDataFile);
else
    fprintf('Error: Could not find the data file at: %s\n', sessionDataFile);
end
```

## Pipeline Parameters: `config/pipeline_config.m`

Global parameters that apply to the entire preprocessing pipeline are defined in the `config/pipeline_config.m` file. This includes key values such as:

*   **`samplingRate`**: The sampling rate of the neural data (e.g., 30000 Hz).
*   **`n_channels_in_dat`**: The number of channels in the processed `.dat` file.

Refer to this file to understand the basic parameters used during data processing.

## Filtering for High-Quality Neurons

The `session_data.spikes` structure contains all the spike data from Kilosort. To ensure you are analyzing high-quality, well-isolated single neurons, you must filter the units based on the information in the `spikes.cluster_info` table.

The most important column for this purpose is `group`. This column contains the quality label assigned during the manual curation step in Phy. You should typically filter for units where `group` is equal to `'good'`. Other labels, such as `'mua'` (multi-unit activity) or `'noise'`, should usually be excluded from single-neuron analysis.

Example (in MATLAB):
```matlab
% Load your session_data.mat file
load('Feynman_08_12_2025_SNc_session_data.mat');

% Find the indices of 'good' clusters
good_cluster_indices = find(strcmp(session_data.spikes.cluster_info.group, 'good'));

% Get the cluster IDs for the 'good' clusters
good_cluster_ids = session_data.spikes.cluster_info.cluster_id(good_cluster_indices);

% Filter spike times to include only 'good' clusters
good_spikes_mask = ismember(session_data.spikes.clusters, good_cluster_ids);
good_spike_times = session_data.spikes.times(good_spikes_mask);
```

## Linking to Behavioral Task Data

The `session_manifest.csv` file provides the link between a `session_data.mat` file and the specific behavioral task that was performed. The `notes` and `experiment_pc_name` columns can be used to identify the task.

Once you have identified the task, you can find a detailed data dictionary for that task's specific variables in the `/docs/task_data_dictionaries/` directory. These dictionaries explain the meaning of the task-specific columns found in the `session_data.trialInfo` table. This is crucial for understanding the behavioral data associated with the neural recordings.
