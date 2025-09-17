%% plot_aggregated_baseline_comparison.m
%
%   Generates a publication-quality summary figure for a single brain area,
%   showing results from the baseline comparison analysis. A separate figure
%   is generated for each alignment event.
%
% INPUTS:
%   aggregated_data   - A struct containing aggregated data for the brain area.
%   brain_area_name   - A string with the name of the brain area (e.g., 'SC').
%
% Author: Jules
% Date: 2025-09-17
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
        data_path = aggregated_data.baseline_comparison.(event_name).(cond_name);
        time_vector = data_path.time_vector;

        % --- Data Processing ---
        p_values = data_path.sig;
        is_sig = p_values < 0.05;

        post_event_fr = data_path.post_event_fr;
        baseline_fr = data_path.baseline_fr;

        increase = is_sig & (post_event_fr > baseline_fr);
        decrease = is_sig & (post_event_fr < baseline_fr);

        prop_increase = sum(increase, 1) / size(p_values, 1);
        prop_decrease = -sum(decrease, 1) / size(p_values, 1);

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
