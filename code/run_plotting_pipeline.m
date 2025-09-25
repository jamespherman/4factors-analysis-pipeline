% run_plotting_pipeline.m
%
% Description:
%   This script is the main entry point for generating all final,
%   aggregated figures. It loads aggregated data, then loops through each
%   brain area ('SC', 'SNc') to generate a complete set of summary
%   figures for each population.
%
% Author:
%   Jules
%
% Date:
%   2025-09-17

% Start timer and provide user feedback
tic;
disp('Starting the plotting pipeline...');

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
aggFileName = fullfile(project_root, 'data', 'processed', ...
    'aggregated_analysis_data.mat');

% Load analysis plan, which is needed for some plots.
disp('Defining analysis plan...');
[~, analysis_plan] = define_task_conditions();

% Load aggregated data
disp('Loading aggregated analysis data...');
load(aggFileName);
disp('Data loaded successfully.');

%% Main Plotting Loop
% This section iterates through each brain area and calls the relevant
% plotting functions to generate a full set of figures for each area.
brain_areas = {'SC', 'SNc'};
all_data = {aggregated_sc_data, aggregated_snc_data};

for i = 1:length(brain_areas)
    brain_area_name = brain_areas{i};
    aggregated_data = all_data{i};

    fprintf('Generating plots for %s...\n', brain_area_name);

    % Generate aggregated ROC comparison plot
    plot_aggregated_roc_comparison(aggregated_data, brain_area_name);

    % Generate aggregated ANOVA plot
    plot_aggregated_anova(aggregated_data, brain_area_name, analysis_plan);

    % Generate aggregated baseline comparison plot
    plot_aggregated_baseline_comparison(aggregated_data, brain_area_name);

    % Generate aggregated behavior plot
    plot_aggregated_behavior(aggregated_data, brain_area_name, analysis_plan);
end

% End timer and provide user feedback
toc;
disp('Plotting pipeline finished.');
