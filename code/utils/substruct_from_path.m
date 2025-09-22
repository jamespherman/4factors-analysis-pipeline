function S = substruct_from_path(path_to_check)
%CREATE_SUBSTRUCT_FROM_PATH Converts a path string to a substruct array.
%   S = create_substruct_from_path(path_to_check) takes a string,
%   path_to_check, where nested structure fields are separated by the
%   platform-specific file separator (e.g., 'field1/field2/field3').
%
%   It returns S, a struct array that can be used with the subsref
%   function to dynamically access the specified nested field.
%
%   Example:
%       % 1. Create a sample nested struct
%       session_data.analysis.baseline_comparison.fixOff.is_high_reward = 
%       'Success!';
%
%       % 2. Define the path to the desired field
%       path_str = fullfile('analysis', 'baseline_comparison', 'fixOff', 
%       'is_high_reward');
%
%       % 3. Generate the substruct array using this function
%       S = create_substruct_from_path(path_str);
%
%       % 4. Use subsref to retrieve the value
%       retrieved_value = subsref(session_data, S)
%       % retrieved_value will be 'Success!'

    % Argument validation
    if ~ischar(path_to_check) && ~isstring(path_to_check)
        error('Input must be a character vector or a string.');
    end

    % Split the input string into a cell array of field names
    % using the system's file separator ('\' on Windows, '/' on Linux/macOS)
    split_fields = strsplit(path_to_check, filesep);

    % Pre-allocate the struct array for efficiency.
    num_fields = numel(split_fields);
    S(num_fields) = struct('type', [], 'subs', []);

    % Loop through each field name and build the substruct array
    for i = 1:num_fields
        S(i).type = '.'; % Specify dot notation for struct field access
        S(i).subs = split_fields{i}; % The name of the specific field
    end
end