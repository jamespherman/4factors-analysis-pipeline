%% plot_aggregated_baseline_comparison.m
%
%   Generates a publication-quality summary figure for a single brain area,
%   showing the proportion of neurons with significant firing rate
%   increases or decreases compared to a baseline period. The output adheres
%   to the specifications in `docs/plotting_requirements.md`.
%
%   A separate figure is generated for each alignment event.
%
% INPUTS:
%
%   aggregated_data - A struct containing aggregated data for one brain area.
%                     This script requires the `baseline_comparison` field,
%                     structured as:
%
%                     .(event_name).(condition_name) = [1 x nSessions] struct array
%
%                     Each struct in the array must contain:
%                       - .sig: [nNeurons x nTimeBins] matrix of significance
%                               results (+1 for increase, -1 for decrease).
%                       - .time_vector: [1 x nTimeBins] vector of time points.
%
%   brain_area_name - A string with the name of the brain area (e.g., 'SC').
%
% Author: Jules
% Date: 2025-09-24
%

function plot_aggregated_baseline_comparison(aggregated_data, brain_area_name)

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Input Validation
if ~isfield(aggregated_data, 'baseline_comparison')
    warning('plot_aggregated_baseline_comparison:no_data', ...
        'Baseline comparison data is missing from the input struct.');
    return;
end

%% Dynamically Discover and Loop Through Alignment Events
event_names = fieldnames(aggregated_data.baseline_comparison);

for i_event = 1:length(event_names)
    event_name = event_names{i_event};

    % Get conditions for the current event
    condition_names = fieldnames(aggregated_data.baseline_comparison.(event_name));
    n_conditions = length(condition_names);

    if n_conditions == 0; continue; end

    %% Figure and Plotting Setup for the Current Event
    fig = figure('Position', [100, 100, 400 * n_conditions, 400], 'Color', 'w');
    h_axes = gobjects(1, n_conditions);

    colors = richColors();
    plot_color = colors(1,:);
    pos_alpha = 0.6; % Face alpha for increase
    neg_alpha = 0.3; % Face alpha for decrease

    %% Main Plotting Loop for Conditions
    for i_cond = 1:n_conditions
        cond_name = condition_names{i_cond};
        session_data = aggregated_data.baseline_comparison.(event_name).(cond_name);

        % --- Data Aggregation from Struct Array ---
        if isempty(session_data)
            h_axes(1, i_cond) = mySubPlot([1, n_conditions, i_cond]);
            title(h_axes(1, i_cond), sprintf('%s (No Data)', strrep(cond_name, '_', ' ')));
            continue;
        end

        all_sig_matrices = arrayfun(@(s) s.sig, session_data, 'UniformOutput', false);
        sig_matrix = cat(1, all_sig_matrices{:});
        time_vector = session_data(1).time_vector;

        % --- Data Processing ---
        if isempty(sig_matrix)
            prop_increase = zeros(1, length(time_vector));
            prop_decrease = zeros(1, length(time_vector));
        else
            prop_increase = mean(sig_matrix == 1, 1, 'omitnan');
            prop_decrease = -mean(sig_matrix == -1, 1, 'omitnan');
        end

        % --- Plotting ---
        h_axes(1, i_cond) = mySubPlot([1, n_conditions, i_cond]);
        hold on;

        % Plot proportions
        h_inc = barStairsFill(time_vector, zeros(size(prop_increase)), prop_increase);
        set(h_inc(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', pos_alpha);
        delete(h_inc(2:3));

        h_dec = barStairsFill(time_vector, zeros(size(prop_decrease)), prop_decrease);
        set(h_dec(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', neg_alpha);
        delete(h_dec(2:3));

        % --- Formatting ---
        title(h_axes(1, i_cond), strrep(cond_name, '_', ' '), 'Interpreter', 'none');
        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
        xlabel_str = sprintf('Time from %s Onset (s)', strrep(event_name, '_', ' '));
        xlabel(xlabel_str);
    end

    %% Figure Cleanup and Final Touches
    if n_conditions > 1
        set(h_axes(1, 2:end), 'YTickLabel', []);
    end

    ylabel(h_axes(1, 1), 'Proportion of Neurons');

    sgtitle_str = sprintf('Baseline Comparison for %s, aligned to %s', ...
        brain_area_name, strrep(event_name, '_', ' '));
    sgtitle(sgtitle_str, 'Interpreter', 'none');

    allAx = findall(fig, 'Type', 'Axes');
    [~, yLims] = outerLims(allAx);
    set(allAx, 'YLim', yLims, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

    % Save figure
    fig_filename = fullfile(figures_dir, ...
        sprintf('aggregated_baseline_comparison_%s_%s.pdf', ...
        brain_area_name, event_name));
    pdfSave(fig_filename, fig.Position(3:4)/72, fig);
end

end
