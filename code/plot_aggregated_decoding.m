%% plot_aggregated_decoding.m
%
%   Generates a comprehensive summary figure visualizing population decoding
%   results. The output adheres to the specifications in
%   `docs/plotting_requirements.md`.
%
% INPUTS:
%   aggregated_data - A struct containing the aggregated decoding results.
%   brain_area_name - A string with the name of the brain area (e.g., 'SC').
%   analysis_plan   - The analysis plan, used to link generalization tests
%                     to their corresponding standard CV tests.
%   session_ids     - A cell array of session IDs for the given brain area.
%
% Author: Jules
% Date: 2025-09-25 (Refactored to be fully plan-driven)
%

function plot_aggregated_decoding(aggregated_data, brain_area_name, analysis_plan, session_ids)
%% Setup Paths & Style
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Input Validation
if ~isfield(aggregated_data, 'population_decoding')
    warning('plot_aggregated_decoding:no_data', ...
        'Population decoding data is missing from the input struct.');
    return;
end

%% Dynamically Discover Tests from the Analysis Plan
testing_plan = analysis_plan.decoding_plan.testing_plan;
% Filter out the standard cross-validated tests to get the generalization tests
gen_tests_plan = testing_plan(~endsWith({testing_plan.test_name}, '_CV'));
n_tests = length(gen_tests_plan);

if n_tests == 0
    disp('No generalization decoding tests found to plot.');
    return;
end

%% Setup Figure and Subplot Layout
n_cols = 4; % Standardized layout for consistency
n_rows = ceil(n_tests / n_cols);
fig = figure('Position', [100, 100, 400 * n_cols, 350 * n_rows], 'Color', 'w');
h_axes = gobjects(n_tests, 1);

%% Identify Monkeys for Marker Styling
monkey_initials = cellfun(@(x) x(1), session_ids, 'UniformOutput', false);
[unique_monkeys, ~, monkey_idx_map] = unique(monkey_initials);
monkey_markers = {'o', 's', '^', 'd', 'p', 'h'};
if length(unique_monkeys) > length(monkey_markers)
    warning('More monkeys than defined markers; some will share markers.');
    marker_map = @(idx) monkey_markers{mod(idx-1, length(monkey_markers)) + 1};
else
    marker_map = @(idx) monkey_markers{idx};
end

%% Main Plotting Loop
for i_test = 1:n_tests
    plan_item = gen_tests_plan(i_test);
    test_name = plan_item.test_name;
    
    h_axes(i_test) = mySubPlot([n_rows, n_cols, i_test]);
    hold on;

    % --- Find Standard and Generalization Data using the analysis_plan ---
    if ~isfield(aggregated_data.population_decoding, test_name)
        title(h_axes(i_test), {strrep(test_name, '_', ' '), ...
            '(Missing Gen. data)'});
        continue;
    end
    gen_data = aggregated_data.population_decoding.(test_name);

    % Use the plan to find the correct training model tag
    train_model_tag = plan_item.train_model_tag;
    standard_test_name = [train_model_tag, '_CV'];
    
    if ~isfield(aggregated_data.population_decoding, standard_test_name)
        warning('plot_aggregated_decoding:missing_cv', ...
            ['Could not find standard CV test "%s" for ...' ...
            'generalization test "%s". Skipping.'], ...
            standard_test_name, test_name);
        title(h_axes(i_test), {strrep(test_name, '_', ' '), ...
            '(Missing CV data)'});
        continue;
    end
    standard_data = aggregated_data.population_decoding.(standard_test_name);

    % Create maps for quick lookup of session data
    standard_map = containers.Map({standard_data.session_id}, ...
        1:length(standard_data));

    % --- Plot Data for Each Session ---
    for i_session = 1:length(gen_data)
        session_id = gen_data(i_session).session_id;

        if ~isKey(standard_map, session_id) || isempty(session_id)
            % Skip if session_id is empty or no standard counterpart
            continue; 
        end

        std_idx = standard_map(session_id);

        x_val = standard_data(std_idx).accuracy;
        x_ci  = standard_data(std_idx).accuracy_ci;
        y_val = gen_data(i_session).accuracy;
        y_ci  = gen_data(i_session).accuracy_ci;

        monkey_initial = session_id(1);
        monkey_idx = find(strcmp(unique_monkeys, monkey_initial));
        if isempty(monkey_idx), continue; end
        marker_style = marker_map(monkey_idx);

        errorbar(x_val, y_val, y_val - y_ci(1), y_ci(2) - y_val, ...
            x_val - x_ci(1), x_ci(2) - x_val, ...
            marker_style, 'MarkerSize', 6, 'LineWidth', 1, ...
            'CapSize', 0, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', 'k');
    end

    hold on
    plot(0.5*[1 1], [0 1], 'k:', 'LineWidth', 0.5);
    plot([0 1],0.5*[1 1], 'k:', 'LineWidth', 0.5);
    set(gca, 'TickDir', 'Out')
    title(strrep(test_name, '_', ' '), 'Interpreter', 'none');
end

%% Figure Cleanup and Final Touches
try
    valid_axes = h_axes(isgraphics(h_axes));
    if isempty(valid_axes), return; end

    set(valid_axes, 'TickDir', 'Out', 'Color', 'none', 'LineWidth', 1, ...
        'XColor', 'k', 'YColor', 'k');

    legend_handles = gobjects(1, length(unique_monkeys));
    for i = 1:length(unique_monkeys)
        legend_handles(i) = plot(nan, nan, marker_map(i), ...
            'MarkerSize', 8, 'LineStyle', 'none', ...
            'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', 'k');
    end
    legend(legend_handles, unique_monkeys, 'Location', 'best', ...
        'Box', 'off');

    big_ax = axes('Position', [0.05 0.05 0.8 0.9], 'Visible', 'off');
    ylabel(big_ax, 'Generalization Accuracy', 'Visible', 'on', ...
        'FontSize', 12, 'FontWeight', 'bold');
    xlabel(big_ax, 'Standard Accuracy', 'Visible', 'on', ...
        'FontSize', 12, 'FontWeight', 'bold');

    sgtitle_str = sprintf('Decoding Generalization in %s', ...
        brain_area_name);
    sgtitle(sgtitle_str, 'FontSize', 16, 'FontWeight', 'bold');

    figures_dir = fullfile(project_root, 'figures', 'decoding');
    if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
    fig_filename = fullfile(figures_dir, ...
        sprintf('aggregated_decoding_generalization_%s.pdf', ...
        brain_area_name));
    pdfSave(fig_filename, fig.Position(3:4)/72, fig);
catch me
    keyboard
end