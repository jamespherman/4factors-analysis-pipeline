%% plot_aggregated_roc_comparison.m
%
%   Generates a publication-quality summary figure for a single brain area,
%   displaying the proportion of neurons with significant preferences
%   between two conditions, based on a bin-by-bin ROC comparison. The
%   output adheres to the specifications in `docs/plotting_requirements.md`.
%
%   A separate figure is generated for each alignment event.
%
% INPUTS:
%
%   aggregated_data - A struct containing aggregated data for one brain area.
%                     This script requires the `roc_comparison` field, which
%                     is structured as:
%
%                     .(event_name).(comp_name) = [1 x nSessions] struct array
%
%                     Each struct in the array must contain:
%                       - .sig: [nNeurons x nTimeBins] matrix of significance
%                               results (+1 for cond2 pref, -1 for cond1 pref).
%                       - .time_vector: [1 x nTimeBins] vector of time points.
%
%   brain_area_name - A string with the name of the brain area (e.g., 'SC').
%
% Author: Jules
% Date: 2025-09-24
%

function plot_aggregated_roc_comparison(aggregated_data, brain_area_name)

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Input Validation
if ~isfield(aggregated_data, 'roc_comparison')
    warning('plot_aggregated_roc_comparison:no_data', ...
        'ROC comparison data is missing from the input struct.');
    return;
end

%% Dynamically Discover and Loop Through Alignment Events
event_names = fieldnames(aggregated_data.roc_comparison);

for i_event = 1:length(event_names)
    event_name = event_names{i_event};

    % Get comparisons for the current event
    comparison_names = fieldnames(aggregated_data.roc_comparison.(event_name));
    n_comparisons = length(comparison_names);

    if n_comparisons == 0; continue; end

    %% Figure and Plotting Setup for the Current Event
    fig = figure('Position', [100, 100, 400 * n_comparisons, 400], 'Color', 'w');
    h_axes = gobjects(1, n_comparisons);

    colors = richColors;
    posColor = colors(7,:); % Color for cond2 > cond1
    negColor = colors(10,:); % Color for cond1 > cond2

    %% Main Plotting Loop for Comparisons
    for i_comp = 1:n_comparisons
        comp_name = comparison_names{i_comp};
        session_data = aggregated_data.roc_comparison.(event_name).(comp_name);

        % --- Data Aggregation from Struct Array ---
        if isempty(session_data)
            % This case should ideally not be reached if the aggregation script
            % does not create empty fields, but as a safeguard:
            h_axes(1, i_comp) = mySubPlot([1, n_comparisons, i_comp]);
            title(h_axes(1, i_comp), sprintf('%s (No Data)', strrep(comp_name, '_', ' ')));
            continue;
        end

        % Concatenate the .sig matrices from each session into a single matrix.
        all_sig_matrices = arrayfun(@(s) s.sig, session_data, 'UniformOutput', false);
        sig_matrix = cat(1, all_sig_matrices{:});

        % Time vector is consistent across sessions; take from the first.
        time_vector = session_data(1).time_vector;

        % --- Data Processing ---
        if isempty(sig_matrix)
            % If no neurons were found across all sessions, plot zeros.
            prop_cond2_pref = zeros(1, length(time_vector));
            prop_cond1_pref = zeros(1, length(time_vector));
        else
            % Calculate proportion of neurons preferring each condition.
            prop_cond2_pref = mean(sig_matrix == 1, 1, 'omitnan');
            prop_cond1_pref = -mean(sig_matrix == -1, 1, 'omitnan');
        end


        % --- Plotting Data ---
        h_axes(1, i_comp) = mySubPlot([1, n_comparisons, i_comp]);
        hold on;

        h1 = barStairsFill(time_vector, zeros(size(prop_cond2_pref)), prop_cond2_pref);
        delete(h1(2));
        set(h1(1), 'FaceColor', posColor, 'EdgeColor', 'none');
        set(h1(3), 'Color', posColor);

        h2 = barStairsFill(time_vector, zeros(size(prop_cond1_pref)), prop_cond1_pref);
        delete(h2(2));
        set(h2(1), 'FaceColor', negColor, 'EdgeColor', 'none');
        set(h2(3), 'Color', negColor);

        title_str = strrep(comp_name, '_', ' ');
        title(h_axes(1, i_comp), title_str, 'Interpreter', 'none');
        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

        xlabel_str = sprintf('Time from %s Onset (s)', strrep(event_name, '_', ' '));
        xlabel(h_axes(1, i_comp), xlabel_str);
    end

    %% Figure Cleanup and Final Touches
    if n_comparisons > 1
        set(h_axes(:, 2:end), 'YTickLabel', []);
    end

    ylabel(h_axes(1, 1), sprintf('Proportion of Neurons (%s)', brain_area_name));

    allAx = findall(fig, 'Type', 'Axes');
    [~, yLims] = outerLims(allAx);
    set(allAx, 'YLim', yLims, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

    sgtitle_str = sprintf('ROC Comparison for %s, aligned to %s', ...
        brain_area_name, strrep(event_name, '_', ' '));
    sgtitle(sgtitle_str, 'Interpreter', 'none');

    % Save figure
    fig_filename = fullfile(figures_dir, ...
        sprintf('aggregated_roc_comparison_%s_%s.pdf', ...
        brain_area_name, event_name));
    pdfSave(fig_filename, fig.Position(3:4)/72, fig);
end

end
