%% analyze_screening_status.m
%
% Standalone script to analyze current screening statistics and identify
% potential issues. Outputs a text report to the console and saves
% summary data.
%
% Author: Claude Code
% Date: 2026-02-04

%% Setup
clear; clc;
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

%% Load Manifest
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);
fprintf('Loaded manifest with %d sessions\n', height(manifest));

%% Aggregate Screening Data
aggregate_data = aggregate_screening_data(manifest, project_root);

%% Generate Report
fprintf('\n');
fprintf('================================================================================\n');
fprintf('              NEURON SCREENING STATUS REPORT\n');
fprintf('              Generated: %s\n', datestr(now));
fprintf('================================================================================\n\n');

%% 1. Overall Summary
fprintf('1. OVERALL SUMMARY\n');
fprintf('--------------------------------------------------------------------------------\n');
n_total = length(aggregate_data.mean_firing_rates);
n_included = sum(aggregate_data.is_included);
n_sc = sum(strcmp(aggregate_data.brain_areas, 'SC'));
n_snc = sum(strcmp(aggregate_data.brain_areas, 'SNc'));
n_sc_included = sum(aggregate_data.is_included & strcmp(aggregate_data.brain_areas, 'SC'));
n_snc_included = sum(aggregate_data.is_included & strcmp(aggregate_data.brain_areas, 'SNc'));

fprintf('Total neurons:           %d\n', n_total);
fprintf('Total included:          %d (%.1f%%)\n', n_included, 100*n_included/n_total);
fprintf('\n');
fprintf('SC neurons:              %d\n', n_sc);
fprintf('SC included:             %d (%.1f%%)\n', n_sc_included, 100*n_sc_included/n_sc);
fprintf('\n');
fprintf('SNc neurons:             %d\n', n_snc);
fprintf('SNc included:            %d (%.1f%%)\n', n_snc_included, 100*n_snc_included/n_snc);
fprintf('\n');

%% 2. Per-Session Summary
fprintf('2. PER-SESSION SUMMARY\n');
fprintf('--------------------------------------------------------------------------------\n');
session_summary = aggregate_data.session_summary;
disp(session_summary(:, {'session_id', 'brain_area', 'n_neurons', 'n_included', 'pct_included'}));
fprintf('\n');

%% 3. Classification Breakdown
fprintf('3. CLASSIFICATION BREAKDOWN\n');
fprintf('--------------------------------------------------------------------------------\n');

% SC classifications
sc_mask = strcmp(aggregate_data.brain_areas, 'SC');
sc_classes = aggregate_data.neuron_classes(sc_mask);
[sc_unique_classes, ~, sc_class_idx] = unique(sc_classes);
sc_class_counts = accumarray(sc_class_idx, 1);
fprintf('SC Neuron Classes:\n');
for i = 1:length(sc_unique_classes)
    fprintf('  %-15s: %d (%.1f%%)\n', sc_unique_classes{i}, sc_class_counts(i), ...
        100*sc_class_counts(i)/sum(sc_mask));
end
fprintf('\n');

% SNc classifications
snc_mask = strcmp(aggregate_data.brain_areas, 'SNc');
snc_classes = aggregate_data.neuron_classes(snc_mask);
[snc_unique_classes, ~, snc_class_idx] = unique(snc_classes);
snc_class_counts = accumarray(snc_class_idx, 1);
fprintf('SNc Neuron Classes:\n');
for i = 1:length(snc_unique_classes)
    fprintf('  %-15s: %d (%.1f%%)\n', snc_unique_classes{i}, snc_class_counts(i), ...
        100*snc_class_counts(i)/sum(snc_mask));
end
fprintf('\n');

%% 4. Exclusion Reasons
fprintf('4. EXCLUSION REASONS\n');
fprintf('--------------------------------------------------------------------------------\n');

excluded_mask = ~aggregate_data.is_included;
excluded_reasons = aggregate_data.exclusion_reasons(excluded_mask);
[unique_reasons, ~, reason_idx] = unique(excluded_reasons);
reason_counts = accumarray(reason_idx, 1);
[sorted_counts, sort_idx] = sort(reason_counts, 'descend');
fprintf('Top exclusion reasons (n=%d excluded):\n', sum(excluded_mask));
for i = 1:min(10, length(unique_reasons))
    idx = sort_idx(i);
    reason = unique_reasons{idx};
    if isempty(reason)
        reason = '(unknown)';
    end
    fprintf('  %3d: %s\n', sorted_counts(i), reason);
end
fprintf('\n');

%% 5. Threshold Analysis
fprintf('5. THRESHOLD ANALYSIS\n');
fprintf('--------------------------------------------------------------------------------\n');

% SC threshold analysis
sc_frs = aggregate_data.mean_firing_rates(sc_mask);
sc_sparsity = aggregate_data.sparsity(sc_mask);
sc_fr_threshold = 3.0;  % Current threshold
sc_sparsity_threshold = 0.80;

fprintf('SC Neurons (FR threshold = %.1f sp/s, Sparsity threshold = %.0f%%):\n', ...
    sc_fr_threshold, sc_sparsity_threshold*100);
fprintf('  FR range:       %.2f - %.2f sp/s (median: %.2f)\n', ...
    min(sc_frs), max(sc_frs), median(sc_frs));
fprintf('  Sparsity range: %.1f%% - %.1f%% (median: %.1f%%)\n', ...
    min(sc_sparsity)*100, max(sc_sparsity)*100, median(sc_sparsity)*100);

% Count marginal neurons (within 20% of threshold)
sc_marginal_fr = sum(sc_frs > sc_fr_threshold*0.8 & sc_frs < sc_fr_threshold*1.2);
sc_marginal_sparse = sum(sc_sparsity > sc_sparsity_threshold*0.9 & sc_sparsity < sc_sparsity_threshold*1.1);
fprintf('  Marginal FR (within 20%% of threshold):      %d neurons\n', sc_marginal_fr);
fprintf('  Marginal sparsity (within 10%% of threshold): %d neurons\n', sc_marginal_sparse);
fprintf('\n');

% SNc threshold analysis
snc_frs = aggregate_data.mean_firing_rates(snc_mask);
snc_sparsity = aggregate_data.sparsity(snc_mask);
snc_fr_threshold = 1.0;  % Current threshold

fprintf('SNc Neurons (FR threshold = %.1f sp/s, Sparsity threshold = %.0f%%):\n', ...
    snc_fr_threshold, sc_sparsity_threshold*100);
fprintf('  FR range:       %.2f - %.2f sp/s (median: %.2f)\n', ...
    min(snc_frs), max(snc_frs), median(snc_frs));
fprintf('  Sparsity range: %.1f%% - %.1f%% (median: %.1f%%)\n', ...
    min(snc_sparsity)*100, max(snc_sparsity)*100, median(snc_sparsity)*100);

snc_marginal_fr = sum(snc_frs > snc_fr_threshold*0.8 & snc_frs < snc_fr_threshold*1.2);
snc_marginal_sparse = sum(snc_sparsity > sc_sparsity_threshold*0.9 & snc_sparsity < sc_sparsity_threshold*1.1);
fprintf('  Marginal FR (within 20%% of threshold):      %d neurons\n', snc_marginal_fr);
fprintf('  Marginal sparsity (within 10%% of threshold): %d neurons\n', snc_marginal_sparse);
fprintf('\n');

%% 6. Stability Analysis
fprintf('6. STABILITY ANALYSIS\n');
fprintf('--------------------------------------------------------------------------------\n');

stability_scores = aggregate_data.stability_scores;
has_stability = ~isnan(stability_scores);
n_with_stability = sum(has_stability);

if n_with_stability > 0
    fprintf('Neurons with stability info: %d/%d (%.1f%%)\n', n_with_stability, n_total, ...
        100*n_with_stability/n_total);
    fprintf('Stability score range: %.2f - %.2f (median: %.2f)\n', ...
        min(stability_scores(has_stability)), max(stability_scores(has_stability)), ...
        median(stability_scores(has_stability)));
    fprintf('Presence ratio range:  %.2f - %.2f (median: %.2f)\n', ...
        min(aggregate_data.presence_ratios(has_stability)), ...
        max(aggregate_data.presence_ratios(has_stability)), ...
        median(aggregate_data.presence_ratios(has_stability)));
    fprintf('Stable period fraction range: %.2f - %.2f (median: %.2f)\n', ...
        min(aggregate_data.stable_period_frac(has_stability)), ...
        max(aggregate_data.stable_period_frac(has_stability)), ...
        median(aggregate_data.stable_period_frac(has_stability)));
else
    fprintf('No stability information available.\n');
end
fprintf('\n');

%% 7. Modulation Analysis
fprintf('7. MODULATION ANALYSIS\n');
fprintf('--------------------------------------------------------------------------------\n');

n_modulated = sum(aggregate_data.is_modulated);
n_sc_modulated = sum(aggregate_data.is_modulated & sc_mask);
n_snc_modulated = sum(aggregate_data.is_modulated & snc_mask);

fprintf('Total modulated: %d/%d (%.1f%%)\n', n_modulated, n_total, 100*n_modulated/n_total);
fprintf('SC modulated:    %d/%d (%.1f%%)\n', n_sc_modulated, n_sc, 100*n_sc_modulated/n_sc);
fprintf('SNc modulated:   %d/%d (%.1f%%)\n', n_snc_modulated, n_snc, 100*n_snc_modulated/n_snc);
fprintf('\n');

%% 8. Potential Issues
fprintf('8. POTENTIAL ISSUES\n');
fprintf('--------------------------------------------------------------------------------\n');

issues_found = false;

% Issue 1: Very low inclusion rate
if n_included/n_total < 0.3
    fprintf('WARNING: Overall inclusion rate is very low (%.1f%%). Consider relaxing thresholds.\n', ...
        100*n_included/n_total);
    issues_found = true;
end

% Issue 2: High variability between sessions
if ~isempty(session_summary)
    pct_range = max(session_summary.pct_included) - min(session_summary.pct_included);
    if pct_range > 40
        fprintf('WARNING: Large variability in inclusion rate across sessions (%.1f%% range).\n', pct_range);
        fprintf('  Highest: %.1f%% (%s)\n', max(session_summary.pct_included), ...
            session_summary.session_id{session_summary.pct_included == max(session_summary.pct_included)});
        fprintf('  Lowest:  %.1f%% (%s)\n', min(session_summary.pct_included), ...
            session_summary.session_id{session_summary.pct_included == min(session_summary.pct_included)});
        issues_found = true;
    end
end

% Issue 3: Many marginal neurons
if sc_marginal_fr > 0.1 * n_sc || snc_marginal_fr > 0.1 * n_snc
    fprintf('NOTE: Many neurons are near FR threshold boundary. Small changes would affect yield significantly.\n');
    issues_found = true;
end

% Issue 4: Missing stability info
if n_with_stability < n_total * 0.9
    fprintf('NOTE: %.1f%% of neurons missing stability information.\n', ...
        100*(1 - n_with_stability/n_total));
    issues_found = true;
end

if ~issues_found
    fprintf('No major issues detected.\n');
end
fprintf('\n');

%% 9. Yield at Alternative Thresholds
fprintf('9. YIELD AT ALTERNATIVE THRESHOLDS (SC)\n');
fprintf('--------------------------------------------------------------------------------\n');

fr_thresholds = [1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0];
sparsity_thresholds = [0.60, 0.70, 0.80, 0.90];

fprintf('SC Yield Matrix (rows=FR threshold, cols=Sparsity threshold):\n');
fprintf('FR\\Sparse');
for sp = sparsity_thresholds
    fprintf('\t%.0f%%', sp*100);
end
fprintf('\n');

for fr = fr_thresholds
    fprintf('%.1f', fr);
    for sp = sparsity_thresholds
        pass_fr = sc_frs >= fr;
        pass_sparse = sc_sparsity <= sp;
        yield = 100 * sum(pass_fr & pass_sparse) / n_sc;
        fprintf('\t%.1f%%', yield);
    end
    fprintf('\n');
end
fprintf('\n');

fprintf('SNc Yield Matrix:\n');
fr_thresholds_snc = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
fprintf('FR\\Sparse');
for sp = sparsity_thresholds
    fprintf('\t%.0f%%', sp*100);
end
fprintf('\n');

for fr = fr_thresholds_snc
    fprintf('%.1f', fr);
    for sp = sparsity_thresholds
        pass_fr = snc_frs >= fr;
        pass_sparse = snc_sparsity <= sp;
        yield = 100 * sum(pass_fr & pass_sparse) / n_snc;
        fprintf('\t%.1f%%', yield);
    end
    fprintf('\n');
end
fprintf('\n');

%% Save Report
fprintf('================================================================================\n');
fprintf('Report complete.\n');
fprintf('================================================================================\n');

% Save aggregate data for further analysis
output_path = fullfile(project_root, 'data', 'screening_analysis_data.mat');
save(output_path, 'aggregate_data');
fprintf('Aggregate data saved to: %s\n', output_path);
