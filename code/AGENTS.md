# Agent Instructions

When creating new MATLAB scripts in this directory, please ensure you add the `utils` directory to the MATLAB path. This is necessary for helper functions to be found.

You can use the following code snippet at the beginning of your script to achieve this:

```matlab
%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
```
This will ensure that the path is added correctly, regardless of the current working directory.
