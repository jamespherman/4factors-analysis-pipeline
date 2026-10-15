function plot_screening_metrics_session(session_data, screening_info, output_dir, session_name)
% PLOT_SCREENING_METRICS_SESSION Generate diagnostic plots for neuron screening metrics.
%
%   PLOT_SCREENING_METRICS_SESSION(session_data, screening_info, output_dir, session_name)
%
%   Creates a multi-page PDF showing screening metric distributions for a single session:
%   - Page 1: FR vs Sparsity scatter plot with thresholds
%   - Page 2: Stability metrics panel
%   - Page 3: Modulation summary
%   - Page 4: Exclusion breakdown
%
%   INPUTS:
%       session_data   - Standard session data structure with metadata
%       screening_info - Struct array from apply_neuron_screening
%       output_dir     - Directory to save output PDF
%       session_name   - Session identifier for filename
%
%   OUTPUT:
%       Saves PDF to: output_dir/{session_name}_screening_metrics.pdf
%
%   Author: Claude Code
%   Date: 2026-01-23

if nargin < 4 || isempty(session_name)
    if isfield(session_data, 'metadata') && isfield(session_data.metadata, 'unique_id')
        session_name = session_data.metadata.unique_id;
    else
        session_name = 'unknown_session';
    end
end

fprintf('plot_screening_metrics_session: Generating screening metric plots for %s...\n', session_name);

% Ensure output directory exists
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Extract brain area for threshold determination
if isfield(session_data, 'metadata') && isfield(session_data.metadata, 'brain_area')
    brain_area = session_data.metadata.brain_area;
else
    brain_area = 'SC';  % Default
end

% Set thresholds based on brain area
if strcmpi(brain_area, 'SNc')
    fr_threshold = 1.0;
else
    fr_threshold = 3.5;
end
sparsity_threshold = 0.70;

% Extract metrics from screening_info
nNeurons = length(screening_info);
mean_frs = [screening_info.mean_firing_rate]';
sparsity = [screening_info.proportion_empty_bins]';
is_included = [screening_info.included]';

% Check for stability info
has_stability = isfield(session_data, 'stability_info') && ~isempty(session_data.stability_info);
if has_stability
    stability_info = session_data.stability_info;
    stability_scores = [stability_info.stability_score]';
    presence_ratios = [stability_info.presence_ratio]';
    n_stable_trials = [stability_info.n_stable_trials]';
    is_stable = [stability_info.is_stable]';
    stable_frs = [stability_info.stable_fr]';

    % Calculate stable period fraction
    nTrials = length(session_data.eventTimes.trialBegin);
    stable_period_frac = n_stable_trials / nTrials;
else
    stability_scores = ones(nNeurons, 1);
    presence_ratios = ones(nNeurons, 1);
    stable_period_frac = ones(nNeurons, 1);
    is_stable = true(nNeurons, 1);
    stable_frs = mean_frs;
end

% Extract modulation info
is_modulated = false(nNeurons, 1);
for i = 1:nNeurons
    if isfield(screening_info(i), 'is_modulated')
        is_modulated(i) = screening_info(i).is_modulated;
    end
end

% Extract exclusion reasons
exclusion_reasons = cell(nNeurons, 1);
for i = 1:nNeurons
    if isfield(screening_info(i), 'exclusion_reason') && ~isempty(screening_info(i).exclusion_reason)
        exclusion_reasons{i} = screening_info(i).exclusion_reason;
    else
        exclusion_reasons{i} = '';
    end
end

% Output filename
output_file = fullfile(output_dir, [session_name, '_screening_metrics.pdf']);

%% Page 1: FR vs Sparsity Scatter Plot
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 800], 'Visible', 'off');

% Main scatter plot
ax = axes('Position', [0.12 0.15 0.75 0.72]);

% Color by inclusion status
colors = zeros(nNeurons, 3);
colors(is_included, :) = repmat([0.2 0.7 0.2], sum(is_included), 1);  % Green for included
colors(~is_included, :) = repmat([0.7 0.2 0.2], sum(~is_included), 1);  % Red for excluded

% Size by stability score
if has_stability
    sizes = 20 + 80 * stability_scores;  % Range 20-100
else
    sizes = 50 * ones(nNeurons, 1);
end

% Plot with log scale for FR
scatter(ax, mean_frs, sparsity, sizes, colors, 'filled', 'MarkerFaceAlpha', 0.7);
set(ax, 'XScale', 'log');
hold on;

% Add threshold lines
xlims = [min(mean_frs(mean_frs > 0)) * 0.5, max(mean_frs) * 2];
ylims = [0, 1];
xlim(xlims);
ylim(ylims);

% FR threshold line
plot([fr_threshold, fr_threshold], ylims, 'b--', 'LineWidth', 2);
text(fr_threshold * 1.1, 0.02, sprintf('FR = %.1f sp/s', fr_threshold), ...
    'Color', 'b', 'FontSize', 10, 'VerticalAlignment', 'bottom');

% Sparsity threshold line
plot(xlims, [sparsity_threshold, sparsity_threshold], 'b--', 'LineWidth', 2);
text(xlims(1) * 1.5, sparsity_threshold + 0.02, sprintf('Sparsity = %.0f%%', sparsity_threshold*100), ...
    'Color', 'b', 'FontSize', 10);

% Count neurons in each quadrant
n_pass_both = sum(mean_frs >= fr_threshold & sparsity <= sparsity_threshold);
n_fail_fr = sum(mean_frs < fr_threshold & sparsity <= sparsity_threshold);
n_fail_sparsity = sum(mean_frs >= fr_threshold & sparsity > sparsity_threshold);
n_fail_both = sum(mean_frs < fr_threshold & sparsity > sparsity_threshold);

% Add quadrant annotations
text(xlims(2) * 0.7, 0.35, sprintf('Pass both: %d', n_pass_both), 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.2 0.7 0.2]);
text(xlims(1) * 2, 0.35, sprintf('Low FR: %d', n_fail_fr), 'FontSize', 11, 'Color', [0.5 0.5 0.5]);
text(xlims(2) * 0.7, 0.85, sprintf('Sparse: %d', n_fail_sparsity), 'FontSize', 11, 'Color', [0.5 0.5 0.5]);
text(xlims(1) * 2, 0.85, sprintf('Both fail: %d', n_fail_both), 'FontSize', 11, 'Color', [0.7 0.2 0.2]);

xlabel('Mean Firing Rate (sp/s)', 'FontSize', 12);
ylabel('Proportion Empty Bins', 'FontSize', 12);
title(sprintf('%s: FR vs Sparsity (N=%d, %d included)', session_name, nNeurons, sum(is_included)), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Add legend
h1 = scatter(nan, nan, 50, [0.2 0.7 0.2], 'filled');
h2 = scatter(nan, nan, 50, [0.7 0.2 0.2], 'filled');
legend([h1, h2], {'Included', 'Excluded'}, 'Location', 'northeast');

% Add marginal histogram for FR (top)
ax_top = axes('Position', [0.12 0.87 0.75 0.08]);
histogram(log10(mean_frs(mean_frs > 0)), 30, 'FaceColor', [0.5 0.5 0.8], 'EdgeColor', 'none');
set(ax_top, 'XLim', log10(xlims), 'XTick', [], 'YTick', []);
box off;

% Add marginal histogram for sparsity (right)
ax_right = axes('Position', [0.87 0.15 0.08 0.72]);
histogram(sparsity, 30, 'Orientation', 'horizontal', 'FaceColor', [0.5 0.5 0.8], 'EdgeColor', 'none');
set(ax_right, 'YLim', ylims, 'XTick', [], 'YTick', []);
box off;

% Save page 1
exportgraphics(fig, output_file, 'ContentType', 'vector');
close(fig);

%% Page 2: Stability Metrics Panel
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 800], 'Visible', 'off');

% 2a. Histogram: Stable period duration
subplot(2, 2, 1);
histogram(stable_period_frac * 100, 20, 'FaceColor', [0.4 0.7 0.9], 'EdgeColor', 'w');
xlabel('Stable Period (% of trials)');
ylabel('Count');
title('Stable Period Duration');
xline(100, 'r--', 'LineWidth', 1.5);
text(95, max(ylim) * 0.9, 'Full session', 'HorizontalAlignment', 'right', 'Color', 'r');

% 2b. Histogram: Presence ratio (included vs excluded)
subplot(2, 2, 2);
hold on;
edges = linspace(0, 1, 21);
histogram(presence_ratios(is_included), edges, 'FaceColor', [0.2 0.7 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(presence_ratios(~is_included), edges, 'FaceColor', [0.7 0.2 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('Presence Ratio');
ylabel('Count');
title('Presence Ratio Distribution');
legend({'Included', 'Excluded'}, 'Location', 'northwest');

% 2c. Histogram: Stability score
subplot(2, 2, 3);
hold on;
edges = linspace(0, 1, 21);
histogram(stability_scores(is_included), edges, 'FaceColor', [0.2 0.7 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
histogram(stability_scores(~is_included), edges, 'FaceColor', [0.7 0.2 0.2], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
xlabel('Stability Score');
ylabel('Count');
title('Stability Score Distribution');
legend({'Included', 'Excluded'}, 'Location', 'northwest');

% 2d. Scatter: Stable FR vs Total FR
subplot(2, 2, 4);
hold on;
scatter(mean_frs(is_included), stable_frs(is_included), 30, [0.2 0.7 0.2], 'filled', 'MarkerFaceAlpha', 0.6);
scatter(mean_frs(~is_included), stable_frs(~is_included), 30, [0.7 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.6);
max_fr = max([mean_frs; stable_frs]) * 1.1;
plot([0 max_fr], [0 max_fr], 'k--', 'LineWidth', 1);
xlabel('Active Period FR (sp/s)');
ylabel('Stable Period FR (sp/s)');
title('FR: Stable vs Active Period');
axis equal;
xlim([0 max_fr]);
ylim([0 max_fr]);
legend({'Included', 'Excluded', 'y=x'}, 'Location', 'southeast');

sgtitle(sprintf('%s: Stability Metrics', session_name), 'FontWeight', 'bold', 'FontSize', 14);

exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
close(fig);

%% Page 3: Modulation Summary
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 800], 'Visible', 'off');

% 3a. Pie chart: Modulation status
subplot(2, 2, 1);
n_mod_included = sum(is_modulated & is_included);
n_mod_excluded = sum(is_modulated & ~is_included);
n_not_mod = sum(~is_modulated);

if n_mod_included + n_mod_excluded + n_not_mod > 0
    pie_data = [n_mod_included, n_mod_excluded, n_not_mod];
    pie_labels = {sprintf('Mod+Incl (%d)', n_mod_included), ...
                  sprintf('Mod+Excl (%d)', n_mod_excluded), ...
                  sprintf('No Mod (%d)', n_not_mod)};
    % Remove zero segments
    valid = pie_data > 0;
    pie(pie_data(valid), pie_labels(valid));
    title('Modulation Status');
    colormap(gca, [0.2 0.7 0.2; 0.9 0.6 0.2; 0.7 0.2 0.2]);
end

% 3b. Bar chart: Modulation by event type (if available)
subplot(2, 2, 2);
% Count modulation events from modulated_at field
event_counts = struct('targetOn', 0, 'saccadeOnset', 0, 'other', 0);
for i = 1:nNeurons
    if isfield(screening_info(i), 'modulated_at') && ~isempty(screening_info(i).modulated_at)
        mod_str = screening_info(i).modulated_at;
        if iscell(mod_str)
            mod_str = strjoin(mod_str, ' ');
        end
        if contains(mod_str, 'targetOn')
            event_counts.targetOn = event_counts.targetOn + 1;
        end
        if contains(mod_str, 'saccadeOnset')
            event_counts.saccadeOnset = event_counts.saccadeOnset + 1;
        end
    end
end
bar([event_counts.targetOn, event_counts.saccadeOnset], 'FaceColor', [0.4 0.6 0.8]);
set(gca, 'XTickLabel', {'Target On', 'Saccade Onset'});
ylabel('Neuron Count');
title('Modulation by Event');

% 3c. Histogram: Number of modulated locations per neuron
subplot(2, 2, 3);
n_mod_locations = zeros(nNeurons, 1);
for i = 1:nNeurons
    if isfield(screening_info(i), 'modulated_at') && ~isempty(screening_info(i).modulated_at)
        if iscell(screening_info(i).modulated_at)
            n_mod_locations(i) = length(screening_info(i).modulated_at);
        elseif ischar(screening_info(i).modulated_at)
            n_mod_locations(i) = length(strsplit(screening_info(i).modulated_at, ','));
        end
    end
end
hold on;
edges = -0.5:1:(max(n_mod_locations)+0.5);
histogram(n_mod_locations(is_included), edges, 'FaceColor', [0.2 0.7 0.2], 'FaceAlpha', 0.6);
histogram(n_mod_locations(~is_included), edges, 'FaceColor', [0.7 0.2 0.2], 'FaceAlpha', 0.6);
xlabel('Number of Modulated Locations');
ylabel('Count');
title('Modulation Location Count');
legend({'Included', 'Excluded'}, 'Location', 'northeast');

% 3d. Text summary panel
subplot(2, 2, 4);
axis off;

passes_fr = sum([screening_info.passes_fr_threshold]);
passes_sparsity = sum([screening_info.passes_sparsity_threshold]);

summary_text = {
    sprintf('Session: %s', session_name), ...
    sprintf('Brain Area: %s', brain_area), ...
    '', ...
    sprintf('Total Neurons: %d', nNeurons), ...
    sprintf('Pass FR (>=%.1f sp/s): %d (%.1f%%)', fr_threshold, passes_fr, 100*passes_fr/nNeurons), ...
    sprintf('Pass Sparsity (<=%.0f%%): %d (%.1f%%)', sparsity_threshold*100, passes_sparsity, 100*passes_sparsity/nNeurons), ...
    sprintf('Show Modulation: %d (%.1f%%)', sum(is_modulated), 100*sum(is_modulated)/nNeurons), ...
    '', ...
    sprintf('Final Included: %d (%.1f%%)', sum(is_included), 100*sum(is_included)/nNeurons)
};

text(0.1, 0.9, summary_text, 'FontSize', 12, 'VerticalAlignment', 'top', ...
    'FontName', 'FixedWidth');

sgtitle(sprintf('%s: Modulation Summary', session_name), 'FontWeight', 'bold', 'FontSize', 14);

exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
close(fig);

%% Page 4: Exclusion Breakdown
fig = figure('Color', 'w', 'Position', [100, 100, 1000, 800], 'Visible', 'off');

% Categorize exclusion reasons
exclusion_categories = {'Low FR', 'Sparse', 'No modulation', 'Low FR & Sparse', ...
                        'Low FR & No mod.', 'No mod. & Sparse', 'Low FR & No mod. & Sparse', 'Other'};
exclusion_counts = zeros(length(exclusion_categories), 1);

for i = 1:nNeurons
    if is_included(i)
        continue;
    end
    reason = exclusion_reasons{i};
    if isempty(reason)
        reason = 'Other';
    end

    % Categorize
    has_low_fr = contains(reason, 'Low FR', 'IgnoreCase', true);
    has_sparse = contains(reason, 'Sparse', 'IgnoreCase', true);
    has_no_mod = contains(reason, 'mod', 'IgnoreCase', true) && contains(reason, 'No', 'IgnoreCase', true);

    if has_low_fr && has_sparse && has_no_mod
        exclusion_counts(7) = exclusion_counts(7) + 1;
    elseif has_low_fr && has_no_mod
        exclusion_counts(5) = exclusion_counts(5) + 1;
    elseif has_no_mod && has_sparse
        exclusion_counts(6) = exclusion_counts(6) + 1;
    elseif has_low_fr && has_sparse
        exclusion_counts(4) = exclusion_counts(4) + 1;
    elseif has_low_fr
        exclusion_counts(1) = exclusion_counts(1) + 1;
    elseif has_sparse
        exclusion_counts(2) = exclusion_counts(2) + 1;
    elseif has_no_mod
        exclusion_counts(3) = exclusion_counts(3) + 1;
    else
        exclusion_counts(8) = exclusion_counts(8) + 1;
    end
end

% 4a. Horizontal bar chart of exclusion reasons
subplot(2, 1, 1);
valid_cats = exclusion_counts > 0;
valid_counts = exclusion_counts(valid_cats);
valid_labels = exclusion_categories(valid_cats);
barh(valid_counts, 'FaceColor', [0.7 0.3 0.3]);
set(gca, 'YTickLabel', valid_labels, 'YDir', 'reverse');
xlabel('Number of Neurons');
title('Exclusion Reasons Breakdown');
xlim([0 max(exclusion_counts) * 1.2 + 1]);

% Add count labels
for k = 1:length(valid_counts)
    text(valid_counts(k) + max(exclusion_counts)*0.02 + 0.5, k, ...
        sprintf('%d', valid_counts(k)), 'VerticalAlignment', 'middle');
end

% 4b. Waterfall/funnel chart showing screening pipeline
subplot(2, 1, 2);

% Calculate cumulative filtering
total = nNeurons;
after_fr = sum([screening_info.passes_fr_threshold]);
after_sparsity = sum([screening_info.passes_fr_threshold] & [screening_info.passes_sparsity_threshold]);
after_mod = sum(is_included);

stages = {'Total', 'Pass FR', 'Pass FR+Sparsity', 'Pass All (Included)'};
counts = [total, after_fr, after_sparsity, after_mod];

bar(counts, 'FaceColor', [0.3 0.6 0.8]);
set(gca, 'XTickLabel', stages);
ylabel('Number of Neurons');
title('Screening Pipeline Funnel');

% Add percentage labels
for j = 1:length(counts)
    text(j, counts(j) + max(counts)*0.03, sprintf('%d (%.0f%%)', counts(j), 100*counts(j)/total), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

% Add arrows showing dropoff
hold on;
for j = 1:length(counts)-1
    dropoff = counts(j) - counts(j+1);
    if dropoff > 0
        text(j + 0.5, (counts(j) + counts(j+1))/2, sprintf('-%d', dropoff), ...
            'HorizontalAlignment', 'center', 'Color', [0.7 0.2 0.2], 'FontSize', 9);
    end
end

sgtitle(sprintf('%s: Exclusion Breakdown', session_name), 'FontWeight', 'bold', 'FontSize', 14);

exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
close(fig);

fprintf('plot_screening_metrics_session: Saved 4 pages to %s\n', output_file);

end
