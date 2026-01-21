%% rerun_roc_comparison.m
%
%   Re-runs only the bin-by-bin ROC comparison analysis for all sessions.
%   Use this when the analysis code has been updated but you don't want
%   to re-run the entire pipeline.
%
%   After running this script, run:
%     1. aggregate_analysis_results.m (to update aggregated data)
%     2. run_plotting_pipeline.m (to regenerate figures)

%% Setup
clear; clc;
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

tic;
giveFeed = @(x) disp([num2str(round(toc, 1)) 's - ' x]);

%% Load Manifest and Analysis Plan
giveFeed('Loading session manifest...');
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

giveFeed('Loading analysis plan...');
[~, analysis_plan] = define_task_conditions();

% Get alignment events
alignment_events = analysis_plan.events;

%% Iterate Through Sessions
for i = 1:height(manifest)
    session_id = manifest.unique_id{i};
    giveFeed(sprintf('--- Processing session: %s ---', session_id));

    % Load session data
    local_processed_dir = fullfile(project_root, 'data', 'processed', session_id);
    local_session_data_path = fullfile(local_processed_dir, ...
        [session_id, '_session_data.mat']);

    if ~exist(local_session_data_path, 'file')
        warning('Session data not found for %s. Skipping.', session_id);
        continue;
    end

    load(local_session_data_path, 'session_data');

    % Ensure metadata is present
    session_data.metadata = table2struct(manifest(i,:));
    session_data.metadata.snc_subregion = parse_snc_subregion(...
        session_data.metadata.grid_hole, session_data.metadata.brain_area);

    % Get conditions and core_data
    conditions = define_task_conditions(session_data);
    core_data = session_data.analysis.core_data;

    % Clear existing ROC comparison results
    if isfield(session_data.analysis, 'roc_comparison')
        session_data.analysis = rmfield(session_data.analysis, 'roc_comparison');
        giveFeed('Cleared existing ROC comparison results.');
    end

    % Re-run ROC comparison for each plan item and event
    for j = 1:length(analysis_plan.roc_plan)
        comp = analysis_plan.roc_plan(j);
        comp_name = comp.name;
        giveFeed(sprintf('--> Running ROC Comparison: %s', comp_name));

        for k = 1:length(alignment_events)
            event_name = alignment_events{k};

            % Create temporary comparison struct with event field
            temp_comp = comp;
            temp_comp.event = event_name;

            % Run the analysis
            result = analyze_roc_comparison(core_data, conditions, temp_comp);

            % Store the result
            session_data.analysis.roc_comparison.(event_name).(comp_name) = result;
        end
    end

    % Save updated session data
    giveFeed('Saving updated session data...');
    save(local_session_data_path, 'session_data', '-v7.3');
    giveFeed(sprintf('--- Finished session: %s ---\n', session_id));
end

giveFeed('ROC comparison re-run complete for all sessions.');
giveFeed('Next steps:');
giveFeed('  1. Run aggregate_analysis_results.m');
giveFeed('  2. Run run_plotting_pipeline.m');
