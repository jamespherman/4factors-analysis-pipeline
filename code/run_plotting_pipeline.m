%% run_plotting_pipeline.m
%
%   Main entry point for generating all aggregated summary figures.
%
%   This script loads the aggregated analysis data and iterates through
%   each brain area ('SC', 'SNc') to generate a complete set of
%   publication-quality figures for each population.
%
%   Generated Figures:
%     - ROC Comparison: Factor selectivity across task epochs
%     - ANOVA: Main effects and interactions from n-way ANOVA
%     - Baseline Comparison: Firing rate comparisons across conditions
%     - Behavior: Behavioral performance summary
%     - Decoding: Generalization analysis across task factors
%     - Window ROC: Epoch-specific selectivity proportions and scatter plots
%     - SNc Subregion: Grouped bar and violin plots comparing rvmSNc vs cdlSNc (SNc only)
%
%   Required Data:
%     - data/processed/aggregated_analysis_data.mat
%
% Author: Jules
% Date: 2025-09-17
% Modified: 2026-01-19 - Added window ROC plotting

%% Initialization
tic;
disp('Starting the plotting pipeline...');

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
aggFileName = fullfile(project_root, 'data', 'processed', ...
    'aggregated_analysis_data.mat');

%% Load Required Data
disp('Defining analysis plan...');
[~, analysis_plan] = define_task_conditions();

disp('Loading aggregated analysis data...');
load(aggFileName);
disp('Data loaded successfully.');

%% Main Plotting Loop
brain_areas = {'SC', 'SNc'};
all_data = {aggregated_sc_data, aggregated_snc_data};

for i = 1:length(brain_areas)
    brain_area_name = brain_areas{i};
    aggregated_data = all_data{i};

    fprintf('Generating plots for %s...\n', brain_area_name);

    % ROC comparison across task epochs
    plot_aggregated_roc_comparison(aggregated_data, brain_area_name);

    % ANOVA main effects and interactions
    plot_aggregated_anova(aggregated_data, brain_area_name, analysis_plan);

    % Baseline firing rate comparisons
    plot_aggregated_baseline_comparison(aggregated_data, brain_area_name);

    % Behavioral performance summary
    plot_aggregated_behavior(aggregated_data, brain_area_name, ...
        analysis_plan);

    % Decoding generalization analysis
    plot_aggregated_decoding(aggregated_data, brain_area_name, ...
        analysis_plan, session_ids.(lower(brain_area_name)));

    % Window-based ROC analysis (proportion significant and scatter plots)
    plot_aggregated_window_roc(aggregated_data, brain_area_name);

    % SNc subregion comparison (only for SNc data)
    if strcmpi(brain_area_name, 'SNc')
        plot_snc_subregion_comparison(aggregated_data, brain_area_name);
    end
end

%% Finish
elapsed_time = toc;
fprintf('Plotting pipeline finished in %.1f seconds.\n', elapsed_time);
