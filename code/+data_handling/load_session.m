function session_data = load_session(unique_id, manifest_path, data_base_path)
% LOAD_SESSION loads a session's data and attaches its metadata.
%
%   session_data = LOAD_SESSION(unique_id, manifest_path, data_base_path)
%   loads the session data for a given unique_id, finds its metadata in
%   the manifest, and returns a single struct with the metadata appended.
%
%   Arguments:
%   - unique_id: A string representing the unique identifier for the session.
%   - manifest_path: A string representing the full path to the manifest CSV file.
%   - data_base_path: A string representing the path to the directory where
%     session data .mat files are stored.
%
%   Returns:
%   - session_data: A struct containing the loaded session data. A new
%     field named 'metadata' is added to this struct, which contains the
%     metadata from the manifest file for the given unique_id.
%
%   Example:
%   session = data_handling.load_session("xyz", "C:\data\manifest.csv", "C:\data\sessions");

% Input validation
if ~ischar(unique_id) && ~isstring(unique_id)
    error('unique_id must be a string.');
end
if ~ischar(manifest_path) && ~isstring(manifest_path)
    error('manifest_path must be a string.');
end
if ~ischar(data_base_path) && ~isstring(data_base_path)
    error('data_base_path must be a string.');
end


% Check if manifest file exists
if ~isfile(manifest_path)
    error('Manifest file not found at: %s', manifest_path);
end

% Read the manifest file
manifest_table = readtable(manifest_path);

% Find the row with the matching unique_id
row_idx = strcmp(manifest_table.unique_id, unique_id);
if ~any(row_idx)
    error('unique_id "%s" not found in the manifest file.', unique_id);
end

% Convert the matching row to a struct
metadata_from_manifest = table2struct(manifest_table(row_idx, :));

% Construct the full path to the data file
data_file_name = [unique_id, '_session_data.mat'];
data_file_path = fullfile(data_base_path, data_file_name);

% Check if data file exists
if ~isfile(data_file_path)
    error('Data file not found at: %s', data_file_path);
end

% Load the session data
loaded_data = load(data_file_path);

% Assuming the struct inside the .mat file is named 'session_data'
if isfield(loaded_data, 'session_data')
    session_data = loaded_data.session_data;
else
    % If the struct has a different name, we might need a more robust way
    % to get it. For now, let's assume we don't know the name and get it
    % dynamically.
    fields = fieldnames(loaded_data);
    if numel(fields) == 1
        session_data = loaded_data.(fields{1});
    else
        error('The .mat file should contain a single struct named session_data.');
    end
end


% Add metadata as a new field
session_data.metadata = metadata_from_manifest;

end
