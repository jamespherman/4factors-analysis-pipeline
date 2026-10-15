function plot_neuron_stability(session_data, stability_info, output_dir, session_name)
% PLOT_NEURON_STABILITY Generate diagnostic plots for neuron stability analysis.
%
%   PLOT_NEURON_STABILITY(session_data, stability_info, output_dir, session_name)
%
%   Creates a multi-page PDF showing stability analysis for each neuron:
%   - Page 1: Summary statistics and distribution plots
%   - Subsequent pages: Per-neuron detail with raster, FR time course, and stability markers
%
%   INPUTS:
%       session_data   - Standard session data structure
%       stability_info - Struct array from detect_stable_trial_period
%       output_dir     - Directory to save output PDF
%       session_name   - (Optional) Session identifier for filename
%
%   OUTPUT:
%       Saves PDF to: output_dir/[session_name]_stability_analysis.pdf
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

fprintf('plot_neuron_stability: Generating stability diagnostic plots for %s...\n', session_name);

% Ensure output directory exists
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Extract data
nNeurons = length(stability_info);
cluster_ids = [stability_info.cluster_id];
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;
trial_start_times = session_data.eventTimes.trialBegin;
trial_end_times = session_data.eventTimes.trialEnd;
nTrials = length(trial_start_times);

% Output filename
output_file = fullfile(output_dir, [session_name, '_stability_analysis.pdf']);

%% Page 1: Summary Statistics
fig = figure('Color', 'w', 'Position', [100, 100, 900, 700], 'Visible', 'off');

% 1a. Distribution of stability scores
subplot(2, 3, 1);
stability_scores = [stability_info.stability_score];
histogram(stability_scores, 20, 'FaceColor', [0.3 0.6 0.9]);
xlabel('Stability Score');
ylabel('Count');
title('Stability Score Distribution');

% 1b. Distribution of presence ratios
subplot(2, 3, 2);
presence_ratios = [stability_info.presence_ratio];
histogram(presence_ratios, 20, 'FaceColor', [0.9 0.5 0.3]);
xlabel('Presence Ratio');
ylabel('Count');
title('Presence Ratio Distribution');

% 1c. Stable vs Unstable pie chart
subplot(2, 3, 3);
n_stable = sum([stability_info.is_stable]);
n_unstable = nNeurons - n_stable;
pie([n_stable, n_unstable], {'Stable', 'Unstable'});
title(sprintf('Stability Classification\n(%d/%d stable)', n_stable, nNeurons));

% 1d. Stable period length distribution
subplot(2, 3, 4);
stable_trials = [stability_info.n_stable_trials];
histogram(stable_trials, 20, 'FaceColor', [0.4 0.8 0.4]);
xlabel('Stable Period (trials)');
ylabel('Count');
title('Stable Period Length');
hold on;
xline(nTrials, 'r--', 'LineWidth', 1.5);
text(nTrials * 0.95, max(ylim) * 0.9, 'Full session', 'HorizontalAlignment', 'right', 'Color', 'r');

% 1e. Stable FR vs Total FR comparison
subplot(2, 3, 5);
stable_frs = [stability_info.stable_fr];
% Calculate total FR for comparison
total_frs = zeros(nNeurons, 1);
total_duration = trial_end_times(end) - trial_start_times(1);
for i = 1:nNeurons
    cluster_spikes = all_spike_times(all_spike_clusters == cluster_ids(i));
    total_frs(i) = length(cluster_spikes) / total_duration;
end
scatter(total_frs, stable_frs, 20, 'filled', 'MarkerFaceAlpha', 0.5);
hold on;
max_fr = max([total_frs; stable_frs(:)]) * 1.1;
plot([0 max_fr], [0 max_fr], 'k--', 'LineWidth', 1);
xlabel('Total Session FR (sp/s)');
ylabel('Stable Period FR (sp/s)');
title('FR: Stable vs Total');
axis equal;
xlim([0 max_fr]);
ylim([0 max_fr]);

% 1f. Text summary
subplot(2, 3, 6);
axis off;
text(0.1, 0.9, sprintf('Session: %s', session_name), 'FontWeight', 'bold', 'FontSize', 12);
text(0.1, 0.75, sprintf('Total neurons: %d', nNeurons));
text(0.1, 0.6, sprintf('Stable neurons: %d (%.1f%%)', n_stable, 100*n_stable/nNeurons));

% Calculate neurons where FR differs significantly
fr_diff = abs(stable_frs(:) - total_frs) ./ max(total_frs, 0.1);
n_diff_10pct = sum(fr_diff > 0.1);
text(0.1, 0.45, sprintf('FR differs >10%%: %d neurons', n_diff_10pct));

mean_stable_length = mean([stability_info.n_stable_trials]);
text(0.1, 0.3, sprintf('Mean stable period: %.1f trials', mean_stable_length));

mean_presence = mean(presence_ratios);
text(0.1, 0.15, sprintf('Mean presence ratio: %.2f', mean_presence));

% Save page 1 using exportgraphics for multi-page PDF support
sgtitle(sprintf('Stability Analysis Summary: %s', session_name), 'FontWeight', 'bold');
exportgraphics(fig, output_file, 'ContentType', 'vector');
close(fig);

%% Per-neuron detail pages (up to 6 neurons per page)
neurons_per_page = 6;
n_pages = ceil(nNeurons / neurons_per_page);

for i_page = 1:n_pages
    fig = figure('Color', 'w', 'Position', [100, 100, 1200, 800], 'Visible', 'off');

    start_idx = (i_page - 1) * neurons_per_page + 1;
    end_idx = min(i_page * neurons_per_page, nNeurons);

    for i_subplot = 1:(end_idx - start_idx + 1)
        i_neuron = start_idx + i_subplot - 1;
        cluster_id = cluster_ids(i_neuron);
        stab = stability_info(i_neuron);

        % Get spike times for this neuron
        cluster_spikes = all_spike_times(all_spike_clusters == cluster_id);

        % Left panel: Sliding window FR
        subplot(neurons_per_page, 2, (i_subplot - 1) * 2 + 1);

        if ~isempty(stab.window_frs)
            window_frs = stab.window_frs;
            nWindows = length(window_frs);
            bar(1:nWindows, window_frs, 'FaceColor', [0.5 0.7 0.9], 'EdgeColor', 'none');
            hold on;

            % Highlight stable region
            if stab.is_stable
                % Estimate which windows correspond to stable period
                % This is approximate since CV calculation uses overlapping windows
                stable_start_frac = stab.stable_trial_start / nTrials;
                stable_end_frac = stab.stable_trial_end / nTrials;
                stable_win_start = max(1, round(stable_start_frac * nWindows));
                stable_win_end = min(nWindows, round(stable_end_frac * nWindows));

                bar(stable_win_start:stable_win_end, window_frs(stable_win_start:stable_win_end), ...
                    'FaceColor', [0.3 0.8 0.3], 'EdgeColor', 'none');
            end

            xlabel('Window');
            ylabel('FR (sp/s)');
        else
            text(0.5, 0.5, 'No window data', 'HorizontalAlignment', 'center');
            axis off;
        end

        % Title with status
        if stab.is_stable
            status_str = sprintf('STABLE (%.0f%% trials)', 100 * stab.n_stable_trials / nTrials);
            title_color = [0 0.5 0];
        else
            status_str = 'UNSTABLE';
            title_color = [0.8 0 0];
        end
        title(sprintf('Cluster %d: %s', cluster_id, status_str), 'Color', title_color);

        % Right panel: Mini raster (trial vs time)
        subplot(neurons_per_page, 2, (i_subplot - 1) * 2 + 2);

        if ~isempty(cluster_spikes)
            % Create trial-aligned raster
            for i_trial = 1:5:nTrials  % Sample every 5th trial for speed
                trial_start = trial_start_times(i_trial);
                trial_end = trial_end_times(i_trial);
                trial_spikes = cluster_spikes(cluster_spikes >= trial_start & cluster_spikes < trial_end);
                trial_spikes_relative = trial_spikes - trial_start;

                if ~isempty(trial_spikes_relative)
                    plot(trial_spikes_relative, i_trial * ones(size(trial_spikes_relative)), ...
                        'k.', 'MarkerSize', 1);
                    hold on;
                end
            end

            % Highlight stable trial range
            if stab.is_stable
                yline(stab.stable_trial_start, 'g-', 'LineWidth', 1.5);
                yline(stab.stable_trial_end, 'g-', 'LineWidth', 1.5);
            end

            xlabel('Time in trial (s)');
            ylabel('Trial');
            ylim([0 nTrials]);
        else
            text(0.5, 0.5, 'No spikes', 'HorizontalAlignment', 'center');
            axis off;
        end

        % Add metrics text
        text_str = sprintf('FR: %.1f sp/s, PR: %.2f', stab.stable_fr, stab.presence_ratio);
        title(text_str, 'FontSize', 8);
    end

    sgtitle(sprintf('Neuron Stability Details - Page %d/%d', i_page, n_pages));

    % Append to PDF using exportgraphics
    exportgraphics(fig, output_file, 'ContentType', 'vector', 'Append', true);
    close(fig);
end

fprintf('plot_neuron_stability: Saved %d pages to %s\n', n_pages + 1, output_file);

end
