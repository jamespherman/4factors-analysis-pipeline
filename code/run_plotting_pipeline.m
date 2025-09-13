% run_plotting_pipeline.m
%
% Description:
%   This script serves as the main entry point for generating all final,
%   aggregated figures for the project. It loads aggregated data and then
%   calls a series of plotting functions to generate the figures.
%
% Author:
%   Jules
%
% Date:
%   2025-09-12

% Start timer and provide user feedback
tic;
disp('Starting the plotting pipeline...');

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
aggFileName = fullfile(project_root, 'data', 'processed', ...
    'aggregated_analysis_data.mat');

% Load aggregated data
disp('Loading aggregated analysis data...');
load(aggFileName);
disp('Data loaded successfully.');

% Call plotting functions
disp('Generating aggregated ROC comparison plot...');
plot_aggregated_roc_comparison(aggregated_sc_data, aggregated_snc_data);

disp('Generating aggregated ANOVA plot...');
plot_aggregated_anova(aggregated_sc_data, aggregated_snc_data);

disp('Generating aggregated baseline comparison plot...');
plot_aggregated_baseline_comparison(aggregated_sc_data, aggregated_snc_data);

% End timer and provide user feedback
toc;
disp('Plotting pipeline finished.');
