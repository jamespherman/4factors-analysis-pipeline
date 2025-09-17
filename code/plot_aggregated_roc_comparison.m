%% plot_aggregated_roc_comparison.m
%
%   Generates a publication-quality summary figure for a single brain area,
%   displaying results from the bin-by-bin ROC comparison analysis. A
%   separate figure is generated for each alignment event.
%
% INPUTS:
%   aggregated_data   - A struct containing aggregated data for one brain area.
%   brain_area_name   - A string with the name of the brain area (e.g., 'SC').
%
% Author: Jules
% Date: 2025-09-17
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
        data_path = aggregated_data.roc_comparison.(event_name).(comp_name);

        % --- Time Vector and Firing Rate Data ---
        time_vector = data_path.time_vector;

        % --- Data Processing ---
        p_values = data_path.sig;
        is_sig = p_values < 0.05;

        fr_cond1 = data_path.cond1_fr;
        fr_cond2 = data_path.cond2_fr;

        pref_cond1 = is_sig & (fr_cond1 > fr_cond2);
        pref_cond2 = is_sig & (fr_cond2 > fr_cond1);

        count_cond1 = -sum(pref_cond1, 1);
        count_cond2 = sum(pref_cond2, 1);

        % --- Plotting Data ---
        h_axes(1, i_comp) = mySubPlot([1, n_comparisons, i_comp]);
        hold on;

        h1 = barStairsFill(time_vector, zeros(size(count_cond2)), count_cond2);
        delete(h1(2));
        set(h1(1), 'FaceColor', posColor, 'EdgeColor', 'none');
        set(h1(3), 'Color', posColor);

        h2 = barStairsFill(time_vector, zeros(size(count_cond1)), count_cond1);
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

    ylabel(h_axes(1, 1), sprintf('Count of Neurons (%s)', brain_area_name));

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
