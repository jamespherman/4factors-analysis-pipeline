%% plot_aggregated_roc_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the bin-by-bin ROC comparison analysis.
%   A separate figure is generated for each alignment event.
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-14
%

function plot_aggregated_roc_comparison(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Input Validation
if ~isfield(aggregated_sc_data, 'roc_comparison') || ...
   ~isfield(aggregated_snc_data, 'roc_comparison')
    warning('plot_aggregated_roc_comparison:no_data', ...
        'ROC comparison data is missing from one or both input structs.');
    return;
end

%% Dynamically Discover and Loop Through Alignment Events
event_names = fieldnames(aggregated_sc_data.roc_comparison);

for i_event = 1:length(event_names)
    event_name = event_names{i_event};

    % Get comparisons for the current event
    comparison_names = fieldnames(aggregated_sc_data.roc_comparison.(event_name));
    n_comparisons = length(comparison_names);

    if n_comparisons == 0; continue; end

    %% Figure and Plotting Setup for the Current Event
    fig = figure('Position', [100, 100, 400 * n_comparisons, 800], 'Color', 'w');
    h_axes = gobjects(2, n_comparisons);

    colors = richColors;
    posColor = colors(7,:); % Color for cond2 > cond1
    negColor = colors(10,:); % Color for cond1 > cond2

    %% Main Plotting Loop for Comparisons
    for i_comp = 1:n_comparisons
        comp_name = comparison_names{i_comp};

        sc_data_path = aggregated_sc_data.roc_comparison.(event_name).(comp_name);
        snc_data_path = aggregated_snc_data.roc_comparison.(event_name).(comp_name);

        % --- Time Vector and Firing Rate Data ---
        time_vector = sc_data_path.time_vector;

        % --- SC Data Processing ---
        p_values_sc = sc_data_path.sig;
        is_sig_sc = p_values_sc < 0.05;

        fr_cond1_sc = sc_data_path.cond1_fr;
        fr_cond2_sc = sc_data_path.cond2_fr;

        pref_cond1_sc = is_sig_sc & (fr_cond1_sc > fr_cond2_sc);
        pref_cond2_sc = is_sig_sc & (fr_cond2_sc > fr_cond1_sc);

        count_sc_cond1 = -sum(pref_cond1_sc, 1);
        count_sc_cond2 = sum(pref_cond2_sc, 1);

        % --- SNc Data Processing ---
        p_values_snc = snc_data_path.sig;
        is_sig_snc = p_values_snc < 0.05;

        fr_cond1_snc = snc_data_path.cond1_fr;
        fr_cond2_snc = snc_data_path.cond2_fr;

        pref_cond1_snc = is_sig_snc & (fr_cond1_snc > fr_cond2_snc);
        pref_cond2_snc = is_sig_snc & (fr_cond2_snc > fr_cond1_snc);

        count_snc_cond1 = -sum(pref_cond1_snc, 1);
        count_snc_cond2 = sum(pref_cond2_snc, 1);

        % --- Plotting SC Data (Top Row) ---
        h_axes(1, i_comp) = mySubPlot([2, n_comparisons, i_comp]);
        hold on;

        h_sc1 = barStairsFill(time_vector, zeros(size(count_sc_cond2)), count_sc_cond2);
        delete(h_sc1(2));
        set(h_sc1(1), 'FaceColor', posColor, 'EdgeColor', 'none');
        set(h_sc1(3), 'Color', posColor);

        h_sc2 = barStairsFill(time_vector, zeros(size(count_sc_cond1)), count_sc_cond1);
        delete(h_sc2(2));
        set(h_sc2(1), 'FaceColor', negColor, 'EdgeColor', 'none');
        set(h_sc2(3), 'Color', negColor);

        title_str = strrep(comp_name, '_', ' ');
        title(h_axes(1, i_comp), title_str, 'Interpreter', 'none');
        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

        % --- Plotting SNc Data (Bottom Row) ---
        h_axes(2, i_comp) = mySubPlot([2, n_comparisons, i_comp + n_comparisons]);
        hold on;

        h_snc1 = barStairsFill(time_vector, zeros(size(count_snc_cond2)), count_snc_cond2);
        delete(h_snc1(2));
        set(h_snc1(1), 'FaceColor', posColor, 'EdgeColor', 'none');
        set(h_snc1(3), 'Color', posColor);

        h_snc2 = barStairsFill(time_vector, zeros(size(count_snc_cond1)), count_snc_cond1);
        delete(h_snc2(2));
        set(h_snc2(1), 'FaceColor', negColor, 'EdgeColor', 'none');
        set(h_snc2(3), 'Color', negColor);

        xlabel_str = sprintf('Time from %s Onset (s)', strrep(event_name, '_', ' '));
        xlabel(h_axes(2, i_comp), xlabel_str);
        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
    end

    %% Figure Cleanup and Final Touches
    set(h_axes(1, :), 'XTickLabel', []);
    if n_comparisons > 1
        set(h_axes(:, 2:end), 'YTickLabel', []);
    end

    ylabel(h_axes(1, 1), 'Count of Neurons (SC)');
    ylabel(h_axes(2, 1), 'Count of Neurons (SNc)');

    allAx = findall(fig, 'Type', 'Axes');
    [~, yLims] = outerLims(allAx);
    set(allAx, 'YLim', yLims, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

    sgtitle(sprintf('ROC Comparison aligned to %s', strrep(event_name, '_', ' ')), ...
        'Interpreter', 'none');

    % Save figure
    fig_filename = fullfile(figures_dir, ...
        sprintf('aggregated_roc_comparison_%s.pdf', event_name));
    pdfSave(fig_filename, fig.Position(3:4)/72, fig);
end

end
