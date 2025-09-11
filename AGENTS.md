# Agent Documentation: Core Project Conventions

This document outlines the core conventions for this analysis project. Adhering to these standards is crucial for ensuring code portability, maintainability, and ease of collaboration.

## 1. Directory Structure: Project vs. Data

A fundamental convention of this project is the strict separation of the **Project Directory** (the Git repository containing code) and the **Data Directory** (containing large `.mat` data files).

### 1.1. Project Directory

*   **Purpose**: Contains all MATLAB code, configuration files, and documentation. This is the directory that is under version control with Git.
*   **How to Find It**: When running any script, you should assume that the current working directory is the root of the Project Directory. You can programmatically get this path using the `pwd` command in MATLAB.

    ```matlab
    % Assumes the script is run from the project root.
    project_root = pwd;
    ```

### 1.2. Data Directory

*   **Purpose**: Stores the large session data files (`_session_data.mat`). This directory is located on OneDrive to facilitate sharing and is **not** part of the Git repository.
*   **How to Find It**: Use the provided utility function `utils.findOneDrive()` to programmatically locate the path to the OneDrive data directory. This function is designed to work across different operating systems and user accounts.

    ```matlab
    % Find the base path for the data directory.
    data_base_path = utils.findOneDrive();
    ```

## 2. Data File Structure

All session data files must follow a standardized structure within the Data Directory. Each session's data is stored in a dedicated sub-folder named with its `unique_id`.

### Example Structure:

```
<Data Directory (from findOneDrive)>
|
|-- 20230101_J01_gSac/
|   |-- 20230101_J01_gSac_session_data.mat
|
|-- 20230102_J02_vProb/
|   |-- 20230102_J02_vProb_session_data.mat
|
|-- ... (other session folders)
```

### Accessing Data

To load a session's data, always use the `data_handling.load_session` function. It correctly constructs the full path to the data file using the conventions described above.

**Correct Usage:**
```matlab
unique_id_to_load = '20230101_J01_gSac';
project_root = pwd;
data_base_path = utils.findOneDrive();
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');

session_data = data_handling.load_session(unique_id_to_load, manifest_path, data_base_path);
```
