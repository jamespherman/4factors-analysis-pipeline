%% aggregate_neuron_metrics.m
%
%   Aggregates per-neuron metrics and PSTH data across all completed
%   sessions, guided by a plan defined in define_metrics_aggregation_plan.m.
%
%   This script loops through sessions marked as 'complete' in the
%   session_manifest.csv, loads their corresponding session_data.mat file,
%   and extracts the specified data. The aggregated results are saved to a

%   single file in the data/processed/ directory.
%
%   Author: Jules
%   Date: 2025-10-01
%

function aggregate_neuron_metrics()
%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Initialization
% Load the plan that defines what metrics to aggregate.
metrics_plan = define_metrics_aggregation_plan();

% Get project root and construct paths to data directories.
project_root = fileparts(script_dir);
data_dir = fullfile(project_root, 'data');
config_dir = fullfile(project_root, 'config');

% Load the session manifest.
manifest_path = fullfile(config_dir, 'session_manifest.csv');
session_manifest = readtable(manifest_path);

%% Initialize Output Structures
% Get column names for the per-neuron metrics tables.
metric_cols = {metrics_plan.per_neuron_metrics.ColumnName};

% Create empty tables for SC and SNc neurons. Use double as a
% default type, which is more efficient than cell.
sc_metrics = table('Size', [0, numel(metric_cols)], ...
    'VariableTypes', repmat({'double'}, 1, numel(metric_cols)), ...
    'VariableNames', metric_cols);
snc_metrics = table('Size', [0, numel(metric_cols)], ...
    'VariableTypes', repmat({'double'}, 1, numel(metric_cols)), ...
    'VariableNames', metric_cols);

% Initialize parent structure for PSTH data.
psth_data.sc = struct();
psth_data.snc = struct();

% Programmatically create fields for each event type.
for i = 1:numel(metrics_plan.psth_aggregation.Events)
    event_name = metrics_plan.psth_aggregation.Events{i};
    psth_data.sc.(event_name) = [];
    psth_data.snc.(event_name) = [];
end

% time vector
time_vectors_captured = false;
psth_time_vectors = struct();
for i = 1:numel(metrics_plan.psth_aggregation.Events)
    event_name = metrics_plan.psth_aggregation.Events{i};
    psth_time_vectors.(event_name) = [];
end

%% Main Aggregation Loop
% Filter manifest for sessions with completed analysis.
complete_sessions = session_manifest(strcmp(...
    session_manifest.analysis_status, 'complete'), :);

% Loop through each completed session.
for i = 1:height(complete_sessions)
    session_id = complete_sessions.unique_id{i};
    fprintf('Processing session: %s\n', session_id);

    % Load the session data file.
    session_data_path = fullfile(data_dir, 'processed', session_id, ...
        sprintf('%s_session_data.mat', session_id));
    if ~exist(session_data_path, 'file')
        warning('Could not find session data for %s. Skipping.', ...
            session_id);
        continue;
    end
    session_data = load(session_data_path);

    % Get the brain area for routing data.
    brain_area = complete_sessions.brain_area{i};

    % Get the actual session_data struct from the loaded file.
    if ~isfield(session_data, 'session_data')
        warning('session_data variable not found in %s. Skipping.', ...
            session_id);
        continue;
    end
    s_data = session_data.session_data;

    % Capture the event-specific time vectors from the first valid session
    if ~time_vectors_captured && isfield(s_data, 'analysis') ...
            && isfield(s_data.analysis, 'core_data')

        for e = 1:numel(metrics_plan.psth_aggregation.Events)
            event = metrics_plan.psth_aggregation.Events{e};
            psth_time_vectors.(event) = ...
                s_data.analysis.core_data.spikes.(event).time_vector;
        end
        time_vectors_captured = true; % Set flag so this doesn't run again
    end

    %% Aggregate Per-Neuron Metrics
    num_metrics = numel(metrics_plan.per_neuron_metrics);
    temp_metric_data = cell(1, num_metrics);

    for m = 1:num_metrics
        metric_info = metrics_plan.per_neuron_metrics(m);
        fields = strsplit(metric_info.SourcePath, '.');
        data_vector = s_data;
        for k = 1:numel(fields)
            data_vector = vertcat(data_vector.(fields{k}));
        end
        temp_metric_data{m} = data_vector(:);
    end

    % Create a temporary table from the collected data columns.
    session_table = table(temp_metric_data{:}, 'VariableNames', ...
        metric_cols);

    % Append to the appropriate master table based on brain area.
    if strcmpi(brain_area, 'SC')
        sc_metrics = [sc_metrics; session_table];
    elseif strcmpi(brain_area, 'SNc')
        snc_metrics = [snc_metrics; session_table];
    end

    %% Aggregate Mean PSTHs
    % Get the logical vector for selecting neurons.
    fields = strsplit(metrics_plan.psth_aggregation.SelectorPath, '.');
    selected_neurons = s_data;
    for k = 1:numel(fields)
        selected_neurons = selected_neurons.(fields{k});
    end

    % Get the base path for the spike data.
    fields = strsplit(metrics_plan.psth_aggregation.SourcePath, '.');
    psth_base_struct = s_data;
    for k = 1:numel(fields)
        psth_base_struct = psth_base_struct.(fields{k});
    end

    % Loop through each event to aggregate PSTHs.
    for e = 1:numel(metrics_plan.psth_aggregation.Events)
        event = metrics_plan.psth_aggregation.Events{e};
        field = metrics_plan.psth_aggregation.DataField;

        if isfield(psth_base_struct, event) && ...
                isfield(psth_base_struct.(event), field)
            rates = psth_base_struct.(event).(field);

            % Filter for selected neurons.
            selected_rates = rates;

            % Get mean PSTH across trials, result is [neurons x bins].
            mean_psth = squeeze(mean(selected_rates, 2, 'omitnan'));

            % Append to the correct output structure.
            if strcmpi(brain_area, 'SC')
                psth_data.sc.(event) = [psth_data.sc.(event); ...
                    mean_psth];
            elseif strcmpi(brain_area, 'SNc')
                psth_data.snc.(event) = [psth_data.snc.(event); ...
                    mean_psth];
            end
        else
            warning('PSTH data for event %s not found in %s.', ...
                event, session_id);
        end
    end
end

%% Save Aggregated Data
% Define the output file path.
output_dir = fullfile(data_dir, 'processed');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
output_path = fullfile(output_dir, 'aggregated_neuron_metrics.mat');

% Save the aggregated data structures to the file.
fprintf('Saving aggregated data to %s\n', output_path);
save(output_path, 'sc_metrics', 'snc_metrics', ...
    'psth_data', 'psth_time_vectors', '-v7.3');
fprintf('Aggregation complete.\n');
end