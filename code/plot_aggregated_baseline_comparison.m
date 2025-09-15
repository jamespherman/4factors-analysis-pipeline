%% plot_aggregated_baseline_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the baseline comparison analysis.
%   A separate figure is generated for each alignment event.
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-14
%

function plot_aggregated_baseline_comparison(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Input Validation
if ~isfield(aggregated_sc_data, 'baseline_comparison') || ...
   ~isfield(aggregated_snc_data, 'baseline_comparison')
    warning('plot_aggregated_baseline_comparison:no_data', ...
        'Baseline comparison data is missing from one or both input structs.');
    return;
end

%% Dynamically Discover and Loop Through Alignment Events
event_names = fieldnames(aggregated_sc_data.baseline_comparison);

for i_event = 1:length(event_names)
    event_name = event_names{i_event};

    % Get conditions for the current event
    condition_names = fieldnames(aggregated_sc_data.baseline_comparison.(event_name));
    n_conditions = length(condition_names);

    if n_conditions == 0; continue; end

    %% Figure and Plotting Setup for the Current Event
    fig = figure('Position', [100, 100, 400 * n_conditions, 600], 'Color', 'w');
    h_axes = gobjects(1, n_conditions);

    colors = richColors();
    sc_color = colors(1,:);
    snc_color = colors(2,:);
    pos_alpha = 0.6; % Face alpha for increase
    neg_alpha = 0.3; % Face alpha for decrease

    %% Main Plotting Loop for Conditions
    for i_cond = 1:n_conditions
        cond_name = condition_names{i_cond};

        sc_data_path = aggregated_sc_data.baseline_comparison.(event_name).(cond_name);
        snc_data_path = aggregated_snc_data.baseline_comparison.(event_name).(cond_name);

        time_vector = sc_data_path.time_vector;

        % --- SC Data Processing ---
        p_values_sc = sc_data_path.sig;
        is_sig_sc = p_values_sc < 0.05;

        post_event_fr_sc = sc_data_path.post_event_fr;
        baseline_fr_sc = sc_data_path.baseline_fr;

        increase_sc = is_sig_sc & (post_event_fr_sc > baseline_fr_sc);
        decrease_sc = is_sig_sc & (post_event_fr_sc < baseline_fr_sc);

        prop_sc_increase = sum(increase_sc, 1) / size(p_values_sc, 1);
        prop_sc_decrease = -sum(decrease_sc, 1) / size(p_values_sc, 1);

        % --- SNc Data Processing ---
        p_values_snc = snc_data_path.sig;
        is_sig_snc = p_values_snc < 0.05;

        post_event_fr_snc = snc_data_path.post_event_fr;
        baseline_fr_snc = snc_data_path.baseline_fr;

        increase_snc = is_sig_snc & (post_event_fr_snc > baseline_fr_snc);
        decrease_snc = is_sig_snc & (post_event_fr_snc < baseline_fr_snc);

        prop_snc_increase = sum(increase_snc, 1) / size(p_values_snc, 1);
        prop_snc_decrease = -sum(decrease_snc, 1) / size(p_values_snc, 1);

        % --- Plotting (SC and SNc on same axes) ---
        h_axes(1, i_cond) = mySubPlot([1, n_conditions, i_cond]);
        hold on;

        % Plot SC proportions
        h_sc_inc = barStairsFill(time_vector, zeros(size(prop_sc_increase)), prop_sc_increase);
        set(h_sc_inc(1), 'FaceColor', sc_color, 'EdgeColor', 'none', 'FaceAlpha', pos_alpha);
        delete(h_sc_inc(2:3));

        h_sc_dec = barStairsFill(time_vector, zeros(size(prop_sc_decrease)), prop_sc_decrease);
        set(h_sc_dec(1), 'FaceColor', sc_color, 'EdgeColor', 'none', 'FaceAlpha', neg_alpha);
        delete(h_sc_dec(2:3));

        % Plot SNc proportions
        h_snc_inc = barStairsFill(time_vector, zeros(size(prop_snc_increase)), prop_snc_increase);
        set(h_snc_inc(1), 'FaceColor', snc_color, 'EdgeColor', 'none', 'FaceAlpha', pos_alpha);
        delete(h_snc_inc(2:3));

        h_snc_dec = barStairsFill(time_vector, zeros(size(prop_snc_decrease)), prop_snc_decrease);
        set(h_snc_dec(1), 'FaceColor', snc_color, 'EdgeColor', 'none', 'FaceAlpha', neg_alpha);
        delete(h_snc_dec(2:3));

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

    sgtitle(sprintf('Baseline Comparison aligned to %s', strrep(event_name, '_', ' ')), ...
        'Interpreter', 'none');

    h_leg_sc = patch(NaN, NaN, sc_color, 'FaceAlpha', pos_alpha);
    h_leg_snc = patch(NaN, NaN, snc_color, 'FaceAlpha', pos_alpha);
    legend([h_leg_sc, h_leg_snc], {'SC', 'SNc'}, 'Location', 'northeast', 'Box', 'off');

    allAx = findall(fig, 'Type', 'Axes');
    set(allAx, 'TickDir', 'Out', 'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

    % Save figure
    fig_filename = fullfile(figures_dir, ...
        sprintf('aggregated_baseline_comparison_%s.pdf', event_name));
    pdfSave(fig_filename, fig.Position(3:4)/72, fig);
end

end
