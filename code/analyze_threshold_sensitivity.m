%% analyze_threshold_sensitivity.m
% Detailed analysis of threshold sensitivity for screening criteria recommendations.
%
% Author: Claude Code
% Date: 2026-01-23

clear; clc; close all;

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

fprintf('=== Threshold Sensitivity Analysis ===\n\n');

%% Load aggregate data
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);
aggregate_data = aggregate_screening_data(manifest, project_root);

%% Separate by brain area
is_SC = strcmpi(aggregate_data.brain_areas, 'SC');
is_SNc = strcmpi(aggregate_data.brain_areas, 'SNc');

sc_frs = aggregate_data.mean_firing_rates(is_SC);
sc_sparsity = aggregate_data.sparsity(is_SC);
sc_modulated = aggregate_data.is_modulated(is_SC);
sc_included = aggregate_data.is_included(is_SC);

snc_frs = aggregate_data.mean_firing_rates(is_SNc);
snc_sparsity = aggregate_data.sparsity(is_SNc);
snc_modulated = aggregate_data.is_modulated(is_SNc);
snc_included = aggregate_data.is_included(is_SNc);

fprintf('=== Current Status ===\n');
fprintf('SC neurons: %d total, %d included (%.1f%%)\n', sum(is_SC), sum(sc_included), 100*mean(sc_included));
fprintf('SNc neurons: %d total, %d included (%.1f%%)\n', sum(is_SNc), sum(snc_included), 100*mean(snc_included));
fprintf('\n');

%% SC Analysis
fprintf('=== SC Neuron Analysis ===\n');
fprintf('Current thresholds: FR >= 3.5 sp/s, Sparsity <= 70%%\n\n');

% FR distribution
fprintf('SC Firing Rate Distribution:\n');
fprintf('  Min: %.2f sp/s\n', min(sc_frs));
fprintf('  25th percentile: %.2f sp/s\n', prctile(sc_frs, 25));
fprintf('  Median: %.2f sp/s\n', median(sc_frs));
fprintf('  75th percentile: %.2f sp/s\n', prctile(sc_frs, 75));
fprintf('  Max: %.2f sp/s\n', max(sc_frs));
fprintf('\n');

% Sparsity distribution
fprintf('SC Sparsity Distribution:\n');
fprintf('  Min: %.1f%%\n', 100*min(sc_sparsity));
fprintf('  25th percentile: %.1f%%\n', 100*prctile(sc_sparsity, 25));
fprintf('  Median: %.1f%%\n', 100*median(sc_sparsity));
fprintf('  75th percentile: %.1f%%\n', 100*prctile(sc_sparsity, 75));
fprintf('  Max: %.1f%%\n', 100*max(sc_sparsity));
fprintf('\n');

% Yield at different FR thresholds
fprintf('SC Yield vs FR Threshold (with sparsity <= 70%%):\n');
fr_thresholds = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];
for fr_t = fr_thresholds
    pass_fr = sc_frs >= fr_t;
    pass_sparsity = sc_sparsity <= 0.70;
    pass_both = pass_fr & pass_sparsity;
    pass_all = pass_both & sc_modulated;
    fprintf('  FR >= %.1f: %d pass FR (%.1f%%), %d pass both (%.1f%%), %d pass all (%.1f%%)\n', ...
        fr_t, sum(pass_fr), 100*mean(pass_fr), sum(pass_both), 100*mean(pass_both), ...
        sum(pass_all), 100*mean(pass_all));
end
fprintf('\n');

% Yield at different sparsity thresholds
fprintf('SC Yield vs Sparsity Threshold (with FR >= 3.5):\n');
sparsity_thresholds = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90];
for sp_t = sparsity_thresholds
    pass_fr = sc_frs >= 3.5;
    pass_sparsity = sc_sparsity <= sp_t;
    pass_both = pass_fr & pass_sparsity;
    pass_all = pass_both & sc_modulated;
    fprintf('  Sparsity <= %.0f%%: %d pass sparsity (%.1f%%), %d pass both (%.1f%%), %d pass all (%.1f%%)\n', ...
        100*sp_t, sum(pass_sparsity), 100*mean(pass_sparsity), sum(pass_both), 100*mean(pass_both), ...
        sum(pass_all), 100*mean(pass_all));
end
fprintf('\n');

% Combined threshold sweep
fprintf('SC Combined Threshold Analysis:\n');
fprintf('FR\\Sparsity |  65%%  |  70%%  |  75%%  |  80%%  |  85%%  |\n');
fprintf('------------|-------|-------|-------|-------|-------|\n');
for fr_t = [2.5, 3.0, 3.5, 4.0]
    fprintf('FR >= %.1f   |', fr_t);
    for sp_t = [0.65, 0.70, 0.75, 0.80, 0.85]
        pass_all = (sc_frs >= fr_t) & (sc_sparsity <= sp_t) & sc_modulated;
        fprintf(' %3d   |', sum(pass_all));
    end
    fprintf('\n');
end
fprintf('\n');

% Marginal neurons analysis
fprintf('SC Marginal Neurons (near current thresholds):\n');
marginal_fr = (sc_frs >= 2.5) & (sc_frs < 3.5);
marginal_sparsity = (sc_sparsity > 0.70) & (sc_sparsity <= 0.80);
fprintf('  FR in [2.5, 3.5): %d neurons\n', sum(marginal_fr));
fprintf('    Of these, modulated: %d (%.1f%%)\n', sum(marginal_fr & sc_modulated), ...
    100*mean(sc_modulated(marginal_fr)));
fprintf('    Of these, sparsity <= 70%%: %d\n', sum(marginal_fr & (sc_sparsity <= 0.70)));
fprintf('    Of these, sparsity <= 80%%: %d\n', sum(marginal_fr & (sc_sparsity <= 0.80)));
fprintf('\n');
fprintf('  Sparsity in (70%%, 80%%]: %d neurons\n', sum(marginal_sparsity));
fprintf('    Of these, modulated: %d (%.1f%%)\n', sum(marginal_sparsity & sc_modulated), ...
    100*mean(sc_modulated(marginal_sparsity)));
fprintf('    Of these, FR >= 3.5: %d\n', sum(marginal_sparsity & (sc_frs >= 3.5)));
fprintf('    Of these, FR >= 3.0: %d\n', sum(marginal_sparsity & (sc_frs >= 3.0)));
fprintf('\n');

% Neurons that would be gained with relaxed thresholds
fprintf('SC Neurons Gained with Relaxed Thresholds:\n');
current = (sc_frs >= 3.5) & (sc_sparsity <= 0.70) & sc_modulated;
proposed_3_80 = (sc_frs >= 3.0) & (sc_sparsity <= 0.80) & sc_modulated;
proposed_2p5_80 = (sc_frs >= 2.5) & (sc_sparsity <= 0.80) & sc_modulated;

fprintf('  Current (FR>=3.5, Sp<=70%%): %d neurons\n', sum(current));
fprintf('  Proposed (FR>=3.0, Sp<=80%%): %d neurons (+%d, +%.1f%%)\n', ...
    sum(proposed_3_80), sum(proposed_3_80) - sum(current), ...
    100*(sum(proposed_3_80) - sum(current))/sum(current));
fprintf('  Proposed (FR>=2.5, Sp<=80%%): %d neurons (+%d, +%.1f%%)\n', ...
    sum(proposed_2p5_80), sum(proposed_2p5_80) - sum(current), ...
    100*(sum(proposed_2p5_80) - sum(current))/sum(current));
fprintf('\n');

%% SNc Analysis
fprintf('=== SNc Neuron Analysis ===\n');
fprintf('Current thresholds: FR >= 1.0 sp/s, Sparsity <= 70%%\n\n');

% FR distribution
fprintf('SNc Firing Rate Distribution:\n');
fprintf('  Min: %.2f sp/s\n', min(snc_frs));
fprintf('  25th percentile: %.2f sp/s\n', prctile(snc_frs, 25));
fprintf('  Median: %.2f sp/s\n', median(snc_frs));
fprintf('  75th percentile: %.2f sp/s\n', prctile(snc_frs, 75));
fprintf('  Max: %.2f sp/s\n', max(snc_frs));
fprintf('\n');

% Sparsity distribution
fprintf('SNc Sparsity Distribution:\n');
fprintf('  Min: %.1f%%\n', 100*min(snc_sparsity));
fprintf('  25th percentile: %.1f%%\n', 100*prctile(snc_sparsity, 25));
fprintf('  Median: %.1f%%\n', 100*median(snc_sparsity));
fprintf('  75th percentile: %.1f%%\n', 100*prctile(snc_sparsity, 75));
fprintf('  Max: %.1f%%\n', 100*max(snc_sparsity));
fprintf('\n');

% SNc yield analysis
fprintf('SNc Yield vs Sparsity Threshold (with FR >= 1.0):\n');
for sp_t = sparsity_thresholds
    pass_fr = snc_frs >= 1.0;
    pass_sparsity = snc_sparsity <= sp_t;
    pass_both = pass_fr & pass_sparsity;
    pass_all = pass_both & snc_modulated;
    fprintf('  Sparsity <= %.0f%%: %d pass both (%.1f%%), %d pass all (%.1f%%)\n', ...
        100*sp_t, sum(pass_both), 100*mean(pass_both), sum(pass_all), 100*mean(pass_all));
end
fprintf('\n');

%% Quality check: Are marginal neurons still good quality?
fprintf('=== Quality Assessment of Marginal SC Neurons ===\n');

% Look at stability of marginal neurons
if isfield(aggregate_data, 'stability_scores')
    stability = aggregate_data.stability_scores(is_SC);

    high_quality = (sc_frs >= 3.5) & (sc_sparsity <= 0.70);
    marginal_quality = ((sc_frs >= 2.5) & (sc_frs < 3.5)) | ((sc_sparsity > 0.70) & (sc_sparsity <= 0.80));
    marginal_quality = marginal_quality & ~high_quality;  % Exclude overlap

    fprintf('Stability scores:\n');
    fprintf('  High quality neurons (current criteria): mean = %.3f, median = %.3f\n', ...
        nanmean(stability(high_quality)), nanmedian(stability(high_quality)));
    fprintf('  Marginal neurons: mean = %.3f, median = %.3f\n', ...
        nanmean(stability(marginal_quality)), nanmedian(stability(marginal_quality)));
    fprintf('\n');
end

%% Final Recommendations
fprintf('=== RECOMMENDATIONS ===\n\n');

% Calculate specific gains
current_sc = sum((sc_frs >= 3.5) & (sc_sparsity <= 0.70) & sc_modulated);
new_sc_3_80 = sum((sc_frs >= 3.0) & (sc_sparsity <= 0.80) & sc_modulated);
new_sc_3_75 = sum((sc_frs >= 3.0) & (sc_sparsity <= 0.75) & sc_modulated);

current_snc = sum((snc_frs >= 1.0) & (snc_sparsity <= 0.70) & snc_modulated);
new_snc_80 = sum((snc_frs >= 1.0) & (snc_sparsity <= 0.80) & snc_modulated);

fprintf('Option 1: Conservative adjustment\n');
fprintf('  SC: FR >= 3.0 sp/s, Sparsity <= 75%%\n');
fprintf('  SNc: FR >= 1.0 sp/s, Sparsity <= 75%%\n');
fprintf('  SC yield: %d -> %d (+%.1f%%)\n', current_sc, new_sc_3_75, 100*(new_sc_3_75-current_sc)/current_sc);
fprintf('\n');

fprintf('Option 2: Moderate adjustment\n');
fprintf('  SC: FR >= 3.0 sp/s, Sparsity <= 80%%\n');
fprintf('  SNc: FR >= 1.0 sp/s, Sparsity <= 80%%\n');
fprintf('  SC yield: %d -> %d (+%.1f%%)\n', current_sc, new_sc_3_80, 100*(new_sc_3_80-current_sc)/current_sc);
fprintf('  SNc yield: %d -> %d (+%.1f%%)\n', current_snc, new_snc_80, 100*(new_snc_80-current_snc)/max(current_snc,1));
fprintf('\n');

% Total neurons
total_current = sum(aggregate_data.is_included);
total_new = new_sc_3_80 + new_snc_80;
fprintf('Total improvement (Option 2):\n');
fprintf('  Current: %d neurons\n', total_current);
fprintf('  Proposed: %d neurons\n', total_new);
fprintf('  Gain: +%d neurons (+%.1f%%)\n', total_new - total_current, 100*(total_new - total_current)/total_current);
