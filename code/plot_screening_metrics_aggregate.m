function plot_screening_metrics_aggregate(aggregate_data, output_dir)
% PLOT_SCREENING_METRICS_AGGREGATE Generate aggregate screening visualizations.
%
%   PLOT_SCREENING_METRICS_AGGREGATE(aggregate_data, output_dir)
%
%   Creates a multi-page PDF showing aggregate screening metrics across all sessions:
%   - Page 1: Screening yield summary by session
%   - Page 2: Population metric distributions
%   - Page 3: Grand scatter plot (FR vs Sparsity, all neurons)
%   - Page 4: Session comparison panel
%   - Page 5: Threshold sensitivity analysis
%   - Page 6: Temporal trends
%
%   INPUTS:
%       aggregate_data - Struct from aggregate_screening_data()
%       output_dir     - Directory to save output PDF
%
%   OUTPUT:
%       Saves PDF to: output_dir/aggregate_screening_metrics.pdf
%
%   Author: Claude Code
%   Date: 2026-01-23

fprintf('plot_screening_metrics_aggregate: Generating aggregate screening visualizations...\n');

% Ensure output directory exists
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

output_file = fullfile(output_dir, 'aggregate_screening_metrics.pdf');

% Extract data
session_summary = aggregate_data.session_summary;
n_sessions = height(session_summary);
nNeurons = length(aggregate_data.mean_firing_rates);

% Identify brain areas
is_SC = strcmpi(aggregate_data.brain_areas, 'SC');
is_SNc = strcmpi(aggregate_data.brain_areas, 'SNc');

%% Page 1: Screening Yield Summary by Session
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 800], 'Visible', 'off');

% Grouped bar chart
subplot(2, 1, 1);

% Prepare data matrix: [n_sessions x 4] for (pass_FR, pass_sparsity, modulated, included)
yield_data = [session_summary.pct_pass_fr, session_summary.pct_pass_sparsity, ...
              session_summary.pct_modulated, session_summary.pct_included];

b = bar(yield_data, 'grouped');
b(1).FaceColor = [0.4 0.6 0.9];
b(2).FaceColor = [0.9 0.6 0.4];
b(3).FaceColor = [0.6 0.4 0.9];
b(4).FaceColor = [0.2 0.7 0.2];

set(gca, 'XTickLabel', session_summary.session_id, 'XTickLabelRotation', 45);
ylabel('Percentage (%)');
title('Screening Metrics by Session', 'FontSize', 14);
legend({'Pass FR', 'Pass Sparsity', 'Modulated', 'Included'}, 'Location', 'eastoutside');
ylim([0 105]);
grid on;

% Summary text
subplot(2, 1, 2);
axis off;

% Create summary table as text
summary_text = sprintf(['AGGREGATE SUMMARY\n' ...
    '================\n' ...
    'Total Sessions: %d\n' ...
    'Total Neurons: %d\n' ...
    'Total Included: %d (%.1f%%)\n\n' ...
    'SC Sessions: %d, Neurons: %d, Included: %d (%.1f%%)\n' ...
    'SNc Sessions: %d, Neurons: %d, Included: %d (%.1f%%)\n\n' ...
    'Mean Yield: %.1f%% (range: %.1f%% - %.1f%%)\n' ...
    'Mean FR (all): %.2f sp/s\n' ...
    'Mean Sparsity (all): %.1f%%\n'], ...
    n_sessions, nNeurons, sum(aggregate_data.is_included), 100*mean(aggregate_data.is_included), ...
    sum(strcmpi(session_summary.brain_area, 'SC')), sum(is_SC), sum(aggregate_data.is_included & is_SC), 100*mean(aggregate_data.is_included(is_SC)), ...
    sum(strcmpi(session_summary.brain_area, 'SNc')), sum(is_SNc), sum(aggregate_data.is_included & is_SNc), 100*mean(aggregate_data.is_included(is_SNc)), ...
    mean(session_summary.pct_included), min(session_summary.pct_included), max(session_summary.pct_included), ...
    mean(aggregate_data.mean_firing_rates), 100*mean(aggregate_data.sparsity));

text(0.1, 0.9, summary_text, 'FontSize', 11, 'FontName', 'FixedWidth', 'VerticalAlignment', 'top');

sgtitle('Screening Yield Summary', 'FontWeight', 'bold', 'FontSize', 16);

exportgraphics(fig, output_file, 'ContentType', 'vector');
close(fig);

%% Page 2: Population Metric Distributions
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 800], 'Visible', 'off');

% 2a. Histogram: Mean firing rate (by brain area)
subplot(2, 3, 1);
hold on;
frs_valid = aggregate_data.mean_firing_rates(aggregate_data.mean_firing_rates > 0);
edges = logspace(log10(min(frs_valid)*0.9), log10(max(frs_valid)*1.1), 31);

frs_SC = aggregate_data.mean_firing_rates(is_SC & aggregate_data.mean_firing_rates > 0);
frs_SNc = aggregate_data.mean_firing_rates(is_SNc & aggregate_data.mean_firing_rates > 0);

histogram(frs_SC, edges, 'FaceColor', [0.3 0.6 0.9], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(frs_SNc, edges, 'FaceColor', [0.9 0.5 0.3], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
set(gca, 'XScale', 'log');
xline(3.5, 'b--', 'LineWidth', 2);
xline(1.0, 'r--', 'LineWidth', 2);
xlabel('Mean FR (sp/s)');
ylabel('Count');
title('Firing Rate Distribution');
legend({'SC', 'SNc', 'SC thresh', 'SNc thresh'}, 'Location', 'northeast');

% 2b. Histogram: Sparsity (by brain area)
subplot(2, 3, 2);
hold on;
edges = linspace(0, 1, 31);
histogram(aggregate_data.sparsity(is_SC), edges, 'FaceColor', [0.3 0.6 0.9], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(aggregate_data.sparsity(is_SNc), edges, 'FaceColor', [0.9 0.5 0.3], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xline(0.70, 'k--', 'LineWidth', 2);
xlabel('Proportion Empty Bins');
ylabel('Count');
title('Sparsity Distribution');
legend({'SC', 'SNc', 'Threshold'}, 'Location', 'northwest');

% 2c. Histogram: Stability score (included vs excluded)
subplot(2, 3, 3);
hold on;
valid_stab = ~isnan(aggregate_data.stability_scores);
edges = linspace(0, 1, 21);
histogram(aggregate_data.stability_scores(valid_stab & aggregate_data.is_included), edges, ...
    'FaceColor', [0.2 0.7 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(aggregate_data.stability_scores(valid_stab & ~aggregate_data.is_included), edges, ...
    'FaceColor', [0.7 0.2 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('Stability Score');
ylabel('Count');
title('Stability Score Distribution');
legend({'Included', 'Excluded'}, 'Location', 'northwest');

% 2d. Histogram: Stable period duration
subplot(2, 3, 4);
hold on;
valid_frac = ~isnan(aggregate_data.stable_period_frac);
edges = linspace(0, 1, 21);
histogram(aggregate_data.stable_period_frac(valid_frac & aggregate_data.is_included), edges, ...
    'FaceColor', [0.2 0.7 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(aggregate_data.stable_period_frac(valid_frac & ~aggregate_data.is_included), edges, ...
    'FaceColor', [0.7 0.2 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('Stable Period (fraction of trials)');
ylabel('Count');
title('Stable Period Duration');
legend({'Included', 'Excluded'}, 'Location', 'northwest');

% 2e. Histogram: Presence ratio
subplot(2, 3, 5);
hold on;
valid_pr = ~isnan(aggregate_data.presence_ratios);
histogram(aggregate_data.presence_ratios(valid_pr & aggregate_data.is_included), edges, ...
    'FaceColor', [0.2 0.7 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(aggregate_data.presence_ratios(valid_pr & ~aggregate_data.is_included), edges, ...
    'FaceColor', [0.7 0.2 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('Presence Ratio');
ylabel('Count');
title('Presence Ratio Distribution');
legend({'Included', 'Excluded'}, 'Location', 'northwest');

% 2f. Box plots: FR by session
subplot(2, 3, 6);
session_labels = categorical(aggregate_data.session_ids);
boxchart(session_labels, aggregate_data.mean_firing_rates, 'MarkerStyle', 'none');
set(gca, 'XTickLabelRotation', 45);
ylabel('Mean FR (sp/s)');
title('FR Distribution by Session');
ylim([0 prctile(aggregate_data.mean_firing_rates, 95) * 1.2]);

sgtitle('Population Metric Distributions', 'FontWeight', 'bold', 'FontSize', 16);

exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
close(fig);

%% Page 3: Grand Scatter Plot
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 900], 'Visible', 'off');

% Main scatter plot
ax_main = axes('Position', [0.12 0.12 0.7 0.7]);
hold on;

% Plot excluded first (background)
excl_SC = ~aggregate_data.is_included & is_SC;
excl_SNc = ~aggregate_data.is_included & is_SNc;
incl_SC = aggregate_data.is_included & is_SC;
incl_SNc = aggregate_data.is_included & is_SNc;

scatter(aggregate_data.mean_firing_rates(excl_SC), aggregate_data.sparsity(excl_SC), ...
    20, [0.7 0.7 0.7], 'o', 'filled', 'MarkerFaceAlpha', 0.3);
scatter(aggregate_data.mean_firing_rates(excl_SNc), aggregate_data.sparsity(excl_SNc), ...
    20, [0.7 0.7 0.7], '^', 'filled', 'MarkerFaceAlpha', 0.3);

% Plot included (foreground)
scatter(aggregate_data.mean_firing_rates(incl_SC), aggregate_data.sparsity(incl_SC), ...
    30, [0.2 0.6 0.2], 'o', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(aggregate_data.mean_firing_rates(incl_SNc), aggregate_data.sparsity(incl_SNc), ...
    30, [0.2 0.6 0.2], '^', 'filled', 'MarkerFaceAlpha', 0.6);

set(gca, 'XScale', 'log');
xlim_vals = [min(aggregate_data.mean_firing_rates(aggregate_data.mean_firing_rates > 0)) * 0.5, ...
             max(aggregate_data.mean_firing_rates) * 2];
xlim(xlim_vals);
ylim([0 1]);

% Threshold lines
xline(3.5, 'b-', 'LineWidth', 2, 'Label', 'SC FR=3.5');
xline(1.0, 'r-', 'LineWidth', 2, 'Label', 'SNc FR=1.0');
yline(0.70, 'k--', 'LineWidth', 2, 'Label', 'Sparsity=70%');

xlabel('Mean Firing Rate (sp/s)', 'FontSize', 12);
ylabel('Proportion Empty Bins', 'FontSize', 12);
title(sprintf('All Neurons: FR vs Sparsity (N=%d, %d included)', nNeurons, sum(aggregate_data.is_included)), ...
    'FontSize', 14);

% Legend
h1 = scatter(nan, nan, 30, [0.2 0.6 0.2], 'o', 'filled');
h2 = scatter(nan, nan, 30, [0.2 0.6 0.2], '^', 'filled');
h3 = scatter(nan, nan, 20, [0.7 0.7 0.7], 'o', 'filled');
h4 = scatter(nan, nan, 20, [0.7 0.7 0.7], '^', 'filled');
legend([h1, h2, h3, h4], {'SC Incl', 'SNc Incl', 'SC Excl', 'SNc Excl'}, 'Location', 'northeast');

% Marginal histogram - FR (top)
ax_top = axes('Position', [0.12 0.82 0.7 0.12]);
frs_log = log10(aggregate_data.mean_firing_rates(aggregate_data.mean_firing_rates > 0));
histogram(frs_log, 40, 'FaceColor', [0.5 0.5 0.8], 'EdgeColor', 'none');
set(gca, 'XLim', log10(xlim_vals), 'XTick', [], 'YTick', []);
box off;
title('FR Distribution');

% Marginal histogram - Sparsity (right)
ax_right = axes('Position', [0.82 0.12 0.12 0.7]);
histogram(aggregate_data.sparsity, 40, 'Orientation', 'horizontal', ...
    'FaceColor', [0.5 0.5 0.8], 'EdgeColor', 'none');
set(gca, 'YLim', [0 1], 'XTick', [], 'YTick', []);
box off;
title('Sparsity');

exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
close(fig);

%% Page 4: Session Comparison Panel
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 800], 'Visible', 'off');

% Heatmap of metrics by session
subplot(2, 1, 1);

metrics_matrix = [session_summary.pct_pass_fr, session_summary.pct_pass_sparsity, ...
                  session_summary.pct_modulated, session_summary.pct_stable, ...
                  session_summary.pct_included];

% Handle NaN in pct_stable
metrics_matrix(isnan(metrics_matrix)) = 0;

imagesc(metrics_matrix);
colormap(gca, parula);
colorbar;
caxis([0 100]);

set(gca, 'YTick', 1:n_sessions, 'YTickLabel', session_summary.session_id);
set(gca, 'XTick', 1:5, 'XTickLabel', {'Pass FR', 'Pass Sparsity', 'Modulated', 'Stable', 'Included'});
title('Screening Metrics Heatmap (% passing)', 'FontSize', 14);

% Add text labels on heatmap
for i = 1:n_sessions
    for j = 1:5
        text(j, i, sprintf('%.0f', metrics_matrix(i,j)), ...
            'HorizontalAlignment', 'center', 'Color', 'w', 'FontWeight', 'bold');
    end
end

% Neuron count bar
subplot(2, 1, 2);
bar_data = [session_summary.n_neurons - session_summary.n_included, session_summary.n_included];
b = barh(bar_data, 'stacked');
b(1).FaceColor = [0.7 0.3 0.3];
b(2).FaceColor = [0.2 0.7 0.2];
set(gca, 'YTick', 1:n_sessions, 'YTickLabel', session_summary.session_id);
xlabel('Number of Neurons');
title('Neuron Selection by Session', 'FontSize', 14);
legend({'Excluded', 'Included'}, 'Location', 'southeast');

sgtitle('Session Comparison', 'FontWeight', 'bold', 'FontSize', 16);

exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
close(fig);

%% Page 5: Threshold Sensitivity Analysis
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 800], 'Visible', 'off');

% 5a. Yield vs FR threshold
subplot(2, 2, 1);
fr_thresholds_test = 0.5:0.25:6;
yield_vs_fr_SC = zeros(size(fr_thresholds_test));
yield_vs_fr_SNc = zeros(size(fr_thresholds_test));

for k = 1:length(fr_thresholds_test)
    thresh = fr_thresholds_test(k);
    yield_vs_fr_SC(k) = 100 * mean(aggregate_data.mean_firing_rates(is_SC) >= thresh);
    yield_vs_fr_SNc(k) = 100 * mean(aggregate_data.mean_firing_rates(is_SNc) >= thresh);
end

plot(fr_thresholds_test, yield_vs_fr_SC, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(fr_thresholds_test, yield_vs_fr_SNc, 'r-^', 'LineWidth', 2, 'MarkerSize', 6);
xline(3.5, 'b--', 'LineWidth', 1.5);
xline(1.0, 'r--', 'LineWidth', 1.5);
xlabel('FR Threshold (sp/s)');
ylabel('% Passing');
title('Yield vs FR Threshold');
legend({'SC', 'SNc', 'SC curr', 'SNc curr'}, 'Location', 'northeast');
grid on;

% 5b. Yield vs sparsity threshold
subplot(2, 2, 2);
sparsity_thresholds_test = 0.4:0.05:0.95;
yield_vs_sparsity = zeros(size(sparsity_thresholds_test));

for k = 1:length(sparsity_thresholds_test)
    thresh = sparsity_thresholds_test(k);
    yield_vs_sparsity(k) = 100 * mean(aggregate_data.sparsity <= thresh);
end

plot(sparsity_thresholds_test * 100, yield_vs_sparsity, 'k-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
xline(70, 'k--', 'LineWidth', 1.5);
xlabel('Sparsity Threshold (%)');
ylabel('% Passing');
title('Yield vs Sparsity Threshold');
grid on;

% 5c. Heatmap: Yield vs both thresholds (SC only)
subplot(2, 2, 3);
fr_grid = 1:0.5:6;
sparsity_grid = 0.5:0.05:0.9;
yield_heatmap = zeros(length(sparsity_grid), length(fr_grid));

for i = 1:length(sparsity_grid)
    for j = 1:length(fr_grid)
        yield_heatmap(i, j) = 100 * mean(aggregate_data.mean_firing_rates(is_SC) >= fr_grid(j) & ...
                                         aggregate_data.sparsity(is_SC) <= sparsity_grid(i));
    end
end

imagesc(fr_grid, sparsity_grid * 100, yield_heatmap);
set(gca, 'YDir', 'normal');
colormap(gca, hot);
colorbar;
hold on;
plot(3.5, 70, 'w*', 'MarkerSize', 15, 'LineWidth', 2);
xlabel('FR Threshold (sp/s)');
ylabel('Sparsity Threshold (%)');
title('SC Yield vs Both Thresholds (%)');

% 5d. Summary text
subplot(2, 2, 4);
axis off;

% Calculate marginal neurons
fr_margin_SC = sum(aggregate_data.mean_firing_rates(is_SC) >= 3.5 * 0.9 & ...
                   aggregate_data.mean_firing_rates(is_SC) < 3.5 * 1.1);
fr_margin_SNc = sum(aggregate_data.mean_firing_rates(is_SNc) >= 1.0 * 0.9 & ...
                    aggregate_data.mean_firing_rates(is_SNc) < 1.0 * 1.1);
sparsity_margin = sum(aggregate_data.sparsity >= 0.63 & aggregate_data.sparsity <= 0.77);

summary_text = {
    'THRESHOLD SENSITIVITY', ...
    '=====================', ...
    '', ...
    'Current Thresholds:', ...
    sprintf('  SC FR: 3.5 sp/s'), ...
    sprintf('  SNc FR: 1.0 sp/s'), ...
    sprintf('  Sparsity: 70%%'), ...
    '', ...
    'Marginal Neurons (+/- 10%%):', ...
    sprintf('  Near SC FR threshold: %d', fr_margin_SC), ...
    sprintf('  Near SNc FR threshold: %d', fr_margin_SNc), ...
    sprintf('  Near sparsity threshold: %d', sparsity_margin), ...
    '', ...
    sprintf('If SC FR = 3.0: %.1f%% pass', 100*mean(aggregate_data.mean_firing_rates(is_SC) >= 3.0)), ...
    sprintf('If SC FR = 4.0: %.1f%% pass', 100*mean(aggregate_data.mean_firing_rates(is_SC) >= 4.0)), ...
    sprintf('If Sparsity = 60%%: %.1f%% pass', 100*mean(aggregate_data.sparsity <= 0.60)), ...
    sprintf('If Sparsity = 80%%: %.1f%% pass', 100*mean(aggregate_data.sparsity <= 0.80))
};

text(0.05, 0.95, summary_text, 'FontSize', 10, 'FontName', 'FixedWidth', 'VerticalAlignment', 'top');

sgtitle('Threshold Sensitivity Analysis', 'FontWeight', 'bold', 'FontSize', 16);

exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
close(fig);

%% Page 6: Temporal Trends
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 600], 'Visible', 'off');

% Sort sessions by date (assuming session_id contains date info)
% For now, use session order as proxy for time

% 6a. Metrics over sessions
subplot(2, 1, 1);
hold on;
x = 1:n_sessions;
plot(x, session_summary.mean_fr, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(x, session_summary.mean_sparsity * 100, 'r-^', 'LineWidth', 1.5, 'MarkerSize', 6);
if any(~isnan(session_summary.mean_stability_score))
    plot(x, session_summary.mean_stability_score * 100, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 6);
    legend({'Mean FR (sp/s)', 'Mean Sparsity (%)', 'Mean Stability (%)'}, 'Location', 'best');
else
    legend({'Mean FR (sp/s)', 'Mean Sparsity (%)'}, 'Location', 'best');
end
set(gca, 'XTick', x, 'XTickLabel', session_summary.session_id, 'XTickLabelRotation', 45);
ylabel('Value');
title('Recording Quality Metrics Over Sessions', 'FontSize', 12);
grid on;

% 6b. Yield over sessions
subplot(2, 1, 2);
hold on;

% Separate by brain area
SC_idx = find(strcmpi(session_summary.brain_area, 'SC'));
SNc_idx = find(strcmpi(session_summary.brain_area, 'SNc'));

if ~isempty(SC_idx)
    plot(SC_idx, session_summary.pct_included(SC_idx), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
end
if ~isempty(SNc_idx)
    plot(SNc_idx, session_summary.pct_included(SNc_idx), 'r-^', 'LineWidth', 2, 'MarkerSize', 8);
end

% Add trend line (all sessions)
p = polyfit(x', session_summary.pct_included, 1);
plot(x, polyval(p, x), 'k--', 'LineWidth', 1);

set(gca, 'XTick', x, 'XTickLabel', session_summary.session_id, 'XTickLabelRotation', 45);
ylabel('% Included');
title('Screening Yield Over Sessions', 'FontSize', 12);
ylim([0 100]);
legend({'SC', 'SNc', 'Trend'}, 'Location', 'best');
grid on;

sgtitle('Temporal Trends', 'FontWeight', 'bold', 'FontSize', 16);

exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
close(fig);

fprintf('plot_screening_metrics_aggregate: Saved 6 pages to %s\n', output_file);

end
