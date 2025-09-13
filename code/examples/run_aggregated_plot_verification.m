%% run_aggregated_plot_verification.m
%
% This script verifies the functionality of plot_aggregated_roc_comparison.m
% by generating the necessary data and calling the function.
%
% Author: Jules
% Date: 2025-09-12
%

%% Setup
clear; close all;
fprintf('Setting up paths...\n');
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
% Go up one level to the 'code' directory, then add 'utils'
base_code_dir = fileparts(script_dir);
addpath(base_code_dir);
addpath(fullfile(base_code_dir, 'utils'));

%% Load Data
fprintf('Loading manifest and aggregating results...\n');
manifest_path = fullfile(base_code_dir, '..', 'config', 'session_manifest.csv');
opts = detectImportOptions(manifest_path);
opts.Delimiter = {','};
opts.VariableTypes = {'string','string','string','string','string','string'};
manifest = readtable(manifest_path, opts);

% Aggregate the data
[aggregated_sc_data, aggregated_snc_data] = aggregate_analysis_results(manifest);

%% Run Plotting Function
fprintf('Generating plot...\n');
if isempty(fieldnames(aggregated_sc_data)) || isempty(fieldnames(aggregated_snc_data))
    error('Aggregation failed to produce data for one or both brain areas.');
end

plot_aggregated_roc_comparison(aggregated_sc_data, aggregated_snc_data);

%% Save Figure
fprintf('Saving figure...\n');
figure_path = fullfile(base_code_dir, '..', 'figures', 'aggregated_roc_comparison.png');
% Ensure the directory exists
[fig_dir, ~, ~] = fileparts(figure_path);
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
saveas(gcf, figure_path);
fprintf('Verification complete. Figure saved to %s\n', figure_path);

close(gcf);
