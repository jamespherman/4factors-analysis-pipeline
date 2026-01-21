function [selected_neurons, screening_info, neuron_class] = screen_snc_neurons(session_data, core_data, condition_defs, project_root)
% screen_snc_neurons - Classifies SNc neurons into 'modulated' or 'excluded'.
%
% This function implements a two-tier classification system for SNc neurons:
%
% MODULATED SNc NEURON (all must be true):
%   (a) Significant change in activity (increase or decrease) for at least
%       3 consecutive samples at ANY of: fixation onset, target onset,
%       fixation offset, peri-saccade, or reward delivery. Significance is
%       determined by arrayROC bootstrap CI (95%) excluding 0.5.
%   (b) Mean FR > 1 sp/s (computed over active period)
%   (c) Not sparse (<=70% of 100ms bins empty across active period)
%
% EXCLUDED:
%   - Any of the three criteria above is not met
%
% INPUTS:
%   session_data   - A struct containing session-specific data
%   core_data      - Binned spike data from prepare_core_data
%   condition_defs - Analysis configuration from define_task_conditions
%   project_root   - The root path of the project directory
%
% OUTPUT:
%   selected_neurons - A logical vector (nClusters x 1) where true indicates
%                      a neuron classified as 'modulated'.
%   screening_info   - A struct array (one per neuron) with detailed screening
%                      results including: mean_firing_rate, proportion_empty_bins,
%                      is_modulated, modulated_at, passes_fr_threshold,
%                      passes_sparsity_threshold, neuron_class, exclusion_reason.
%   neuron_class     - A cell array (nClusters x 1) with classification:
%                      'modulated' or 'excluded'.
%
% See also: apply_neuron_screening, compute_task_modulation, screen_sc_neurons

fprintf('screen_snc_neurons: Classifying SNc neurons...\n');

%% Extract configuration
neuron_inclusion = condition_defs.neuron_inclusion;
use_strict_screening = neuron_inclusion.use_strict_screening;
fr_threshold = neuron_inclusion.min_firing_rate_snc;  % 1 sp/s for SNc
max_proportion_empty = neuron_inclusion.max_proportion_empty_bins;  % 0.7 (70%)

%% Get neuron information
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nNeurons = numel(cluster_ids);
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

fprintf('screen_snc_neurons: FR threshold=%.1f sp/s, max empty bins=%.0f%%, %d neurons\n', ...
    fr_threshold, max_proportion_empty * 100, nNeurons);

%% Compute total recording duration from trial times (for reference)
recording_start = session_data.eventTimes.trialBegin(1);
recording_end = session_data.eventTimes.trialEnd(end);
recording_duration = recording_end - recording_start;

%% Initialize outputs
selected_neurons = false(nNeurons, 1);
neuron_class = repmat({'excluded'}, nNeurons, 1);

screening_info = struct();
for i = 1:nNeurons
    screening_info(i).cluster_id = cluster_ids(i);
    screening_info(i).included = false;
    screening_info(i).neuron_class = 'excluded';
    screening_info(i).mean_firing_rate = 0;
    screening_info(i).proportion_empty_bins = 1;
    screening_info(i).passes_fr_threshold = false;
    screening_info(i).passes_sparsity_threshold = false;
    screening_info(i).is_modulated = false;
    screening_info(i).modulated_at = {};
    screening_info(i).modulation_details = struct();
    screening_info(i).exclusion_reason = '';
    screening_info(i).active_start = NaN;
    screening_info(i).active_end = NaN;
    screening_info(i).active_duration = 0;
    screening_info(i).fr_threshold = fr_threshold;
    screening_info(i).sparsity_threshold = max_proportion_empty;
    screening_info(i).recording_duration = recording_duration;
end

%% If use_strict_screening is false, include all neurons
if ~use_strict_screening
    fprintf('screen_snc_neurons: use_strict_screening=false, including all %d neurons.\n', nNeurons);
    selected_neurons = true(nNeurons, 1);
    neuron_class = repmat({'modulated'}, nNeurons, 1);
    for i = 1:nNeurons
        screening_info(i).included = true;
        screening_info(i).neuron_class = 'modulated';
    end
    return;
end

%% Criterion 1 & 3: Mean Firing Rate and Sparsity (computed over ACTIVE PERIOD)
for i = 1:nNeurons
    cluster_id = cluster_ids(i);
    spike_times = all_spike_times(all_spike_clusters == cluster_id);
    n_spikes = length(spike_times);

    if n_spikes < 2
        % Not enough spikes to define an active period
        screening_info(i).mean_firing_rate = 0;
        screening_info(i).proportion_empty_bins = 1;
        screening_info(i).passes_fr_threshold = false;
        screening_info(i).passes_sparsity_threshold = false;
        screening_info(i).active_start = NaN;
        screening_info(i).active_end = NaN;
        screening_info(i).active_duration = 0;
        continue;
    end

    % Define active period for this neuron (first spike to last spike)
    active_start = spike_times(1);
    active_end = spike_times(end);
    active_duration = active_end - active_start;

    % Store active period info
    screening_info(i).active_start = active_start;
    screening_info(i).active_end = active_end;
    screening_info(i).active_duration = active_duration;

    % Handle edge case: very short active period
    if active_duration < 1  % Less than 1 second
        screening_info(i).mean_firing_rate = n_spikes;  % Treat as if 1 second
        screening_info(i).proportion_empty_bins = 0;    % Can't compute meaningfully
        screening_info(i).passes_fr_threshold = n_spikes >= fr_threshold;
        screening_info(i).passes_sparsity_threshold = true;  % Benefit of doubt
        continue;
    end

    % Compute mean firing rate over active period only
    mean_fr = n_spikes / active_duration;
    screening_info(i).mean_firing_rate = mean_fr;
    screening_info(i).passes_fr_threshold = mean_fr >= fr_threshold;

    % Compute sparsity over active period only (100ms bins)
    edges = active_start : 0.1 : active_end;
    if length(edges) < 2
        % Not enough bins
        screening_info(i).proportion_empty_bins = 0;
        screening_info(i).passes_sparsity_threshold = true;
    else
        counts = histcounts(spike_times, edges);
        proportion_empty = sum(counts == 0) / numel(counts);
        screening_info(i).proportion_empty_bins = proportion_empty;
        screening_info(i).passes_sparsity_threshold = proportion_empty <= max_proportion_empty;
    end
end

n_pass_fr = sum([screening_info.passes_fr_threshold]);
n_pass_sparsity = sum([screening_info.passes_sparsity_threshold]);
fprintf('screen_snc_neurons: %d/%d neurons pass FR threshold (>= %.1f sp/s)\n', ...
    n_pass_fr, nNeurons, fr_threshold);
fprintf('screen_snc_neurons: %d/%d neurons pass sparsity threshold (<=%.0f%% empty)\n', ...
    n_pass_sparsity, nNeurons, max_proportion_empty * 100);

%% Criterion 2: Task Modulation (arrayROC-based with 3 consecutive significant bins)
% Call compute_task_modulation which checks all events in modulation_windows:
% fixOn, targetOn, fixOff, saccadeOnset, reward
modulation_info = compute_task_modulation(session_data, core_data, condition_defs);

for i = 1:nNeurons
    screening_info(i).is_modulated = modulation_info(i).is_modulated;
    screening_info(i).modulated_at = modulation_info(i).modulated_at;
    screening_info(i).modulation_details = modulation_info(i).modulation_details;

    % Copy bin_significance if available
    if isfield(modulation_info(i), 'bin_significance')
        screening_info(i).bin_significance = modulation_info(i).bin_significance;
    end
end

n_modulated = sum([modulation_info.is_modulated]);
fprintf('screen_snc_neurons: %d/%d neurons show task modulation (>=3 consecutive sig. bins)\n', ...
    n_modulated, nNeurons);

%% Combined Classification Logic (Three-Part Test)
% Apply all three criteria: FR threshold AND task modulation AND sparsity
for i = 1:nNeurons
    pass_fr = screening_info(i).passes_fr_threshold;
    is_mod = screening_info(i).is_modulated;
    pass_sparse = screening_info(i).passes_sparsity_threshold;

    if pass_fr && is_mod && pass_sparse
        % All three criteria met -> MODULATED
        screening_info(i).included = true;
        screening_info(i).neuron_class = 'modulated';
        screening_info(i).exclusion_reason = '';
        selected_neurons(i) = true;
        neuron_class{i} = 'modulated';
    else
        % Build exclusion reason dynamically
        screening_info(i).included = false;
        screening_info(i).neuron_class = 'excluded';
        failed_criteria = {};

        if ~pass_fr
            failed_criteria{end+1} = sprintf('Low FR (%.2f sp/s)', screening_info(i).mean_firing_rate);
        end
        if ~is_mod
            failed_criteria{end+1} = 'No modulation';
        end
        if ~pass_sparse
            failed_criteria{end+1} = sprintf('Sparse (%.0f%%)', screening_info(i).proportion_empty_bins * 100);
        end

        screening_info(i).exclusion_reason = strjoin(failed_criteria, ' & ');
        selected_neurons(i) = false;
        neuron_class{i} = 'excluded';
    end
end

%% Summary
n_selected = sum(selected_neurons);
fprintf('\n=== SNc Neuron Classification Results ===\n');
fprintf('Modulated:  %d\n', n_selected);
fprintf('Excluded:   %d\n', nNeurons - n_selected);
fprintf('Total:      %d\n', nNeurons);

% Count exclusion reasons
exclusion_reasons = {screening_info.exclusion_reason};
non_empty_reasons = exclusion_reasons(~cellfun(@isempty, exclusion_reasons));
unique_reasons = unique(non_empty_reasons);

if ~isempty(unique_reasons)
    fprintf('\nExclusion breakdown:\n');
    for i = 1:length(unique_reasons)
        reason = unique_reasons{i};
        count = sum(strcmp(exclusion_reasons, reason));
        fprintf('  %s: %d\n', reason, count);
    end
end

%% Generate Diagnostic Figure
if ~isempty(project_root)
    output_dir = fullfile(project_root, 'figures');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fig = figure('Color', 'w', 'Position', [100, 100, 900, 400]);

    % Panel 1: FR vs Sparsity scatter with classification
    subplot(1, 2, 1);
    hold on;

    mean_frs = [screening_info.mean_firing_rate];
    sparsity = [screening_info.proportion_empty_bins] * 100;
    is_modulated = selected_neurons;

    % Plot excluded neurons
    scatter(mean_frs(~is_modulated), sparsity(~is_modulated), ...
        50, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.5);
    % Plot modulated neurons
    scatter(mean_frs(is_modulated), sparsity(is_modulated), ...
        50, [0.2 0.6 0.2], 'filled', 'MarkerFaceAlpha', 0.8);

    % Add threshold lines
    xline(fr_threshold, '--r', 'LineWidth', 1.5);
    yline(max_proportion_empty * 100, '--b', 'LineWidth', 1.5);

    xlabel('Mean Firing Rate (sp/s)');
    ylabel('Sparsity (% empty 100ms bins)');
    title(sprintf('SNc Neuron Screening: %s', session_data.metadata.unique_id), ...
        'Interpreter', 'none');
    legend({'Excluded', 'Modulated', 'FR threshold', 'Sparsity threshold'}, ...
        'Location', 'northeast');
    grid on;
    set(gca, 'XScale', 'log');
    xlim([0.1, max(100, max(mean_frs) * 1.2)]);

    % Panel 2: Modulation event breakdown for selected neurons
    subplot(1, 2, 2);

    % Count modulation by event
    event_names = {'fixOn', 'targetOn', 'fixOff', 'saccadeOnset', 'reward'};
    event_counts = zeros(1, length(event_names));

    for i = 1:nNeurons
        if selected_neurons(i)
            modulated_at = screening_info(i).modulated_at;
            for j = 1:length(modulated_at)
                mod_str = modulated_at{j};
                for k = 1:length(event_names)
                    if contains(mod_str, event_names{k})
                        event_counts(k) = event_counts(k) + 1;
                        break;
                    end
                end
            end
        end
    end

    bar(event_counts, 'FaceColor', [0.2 0.6 0.2]);
    set(gca, 'XTickLabel', {'Fix On', 'Target On', 'Fix Off', 'Saccade', 'Reward'});
    xlabel('Event');
    ylabel('# Neurons Modulated');
    title('Modulation by Event (selected neurons)');
    grid on;

    % Save figure
    figFileName = fullfile(output_dir, [session_data.metadata.unique_id, ...
        '_snc_screening.pdf']);
    try
        pdfSave(figFileName, fig.Position(3:4)/72, fig);
        fprintf('screen_snc_neurons: Diagnostic plot saved to %s\n', figFileName);
    catch ME
        warning('screen_snc_neurons: Could not save figure: %s', ME.message);
    end
    close(fig);
end

fprintf('screen_snc_neurons: Finished. Selected %d/%d neurons as modulated.\n', ...
    n_selected, nNeurons);

end
