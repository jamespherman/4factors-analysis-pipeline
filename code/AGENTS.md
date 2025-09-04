# Agent Instructions

This document provides a set of guidelines and conventions to follow when working with the MATLAB code in this repository. Adhering to these standards will help ensure the code is clean, readable, and maintainable.

## Path Setup

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

## Code Formatting and Style

### Line Length
Each line of code should be limited to a maximum of 75 characters. For longer lines, use an ellipsis (`...`) to wrap to the next line. This improves readability and prevents horizontal scrolling.

### Commenting
It is strongly preferred that comments are added on the line *above* the code to which they refer, rather than to the right of the line.

Good:
```matlab
% This is a descriptive comment about the next line of code.
myVariable = someFunction(parameter1, parameter2);
```

Avoid:
```matlab
myVariable = someFunction(parameter1, parameter2); % This is a less-preferred comment.
```

## Naming Conventions

-   **File Names:** Use `snake_case` for all `.m` file names (e.g., `my_analysis_script.m`).
-   **Function Names:** Use `camelCase` for function names (e.g., `calculateSpikeRate`).
-   **Variable Names:** Use descriptive `camelCase` for variable names (e.g., `spikeCount`, `eventTimes`).

## Function Design

### Function Headers
All functions should include a comprehensive header that serves as its documentation. The header should be accessible via the `help` command in MATLAB and should include:
- A brief, one-line summary of the function's purpose.
- A `Usage` section with an example of how to call the function.
- A detailed description of all `Inputs` and `Outputs`.

Example:
```matlab
function myOutput = myFunction(myInput)
% MYFUNCTION is a brief, one-line summary of what the function does.
%
%   Usage:
%   myOutput = myFunction(myInput)
%
%   Inputs:
%   myInput - A description of this input.
%
%   Outputs:
%   myOutput - A description of what is returned.

... function code ...
end
```

### Input Parsing
For functions with multiple or optional arguments, use the `inputParser` to handle argument validation and default values. This makes the function more robust and easier to use.
