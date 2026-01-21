%% plot_aggregated_decoding.m
%
%   Generates a comprehensive summary figure visualizing population decoding
%   results. The output adheres to the specifications in
%   `docs/plotting_requirements.md`.
%
%   The figure is organized into three column groups using viewport
%   partitioning (mySubPlot Pattern 3):
%     - Visual Epoch: 4x4 grid (training factor × test factor)
%     - Saccade Epoch: 4x4 grid (training factor × test factor)
%     - Cross-Epoch: 4x2 grid (Vis→Sac and Sac→Vis generalization)
%
%   Within each epoch group, columns represent training factors (same X-axis
%   values) and rows represent test factors (generalization target).
%   Diagonal cells (self-tests) and Salience×Identity cells are left empty.
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
% Updated: 2026-01-21 (Reorganized layout with viewport partitioning)
%

function plot_aggregated_decoding(aggregated_data, brain_area_name, analysis_plan, session_ids)
%% Setup Paths & Style
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
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

if isempty(gen_tests_plan)
    disp('No generalization decoding tests found to plot.');
    return;
end

%% Categorize Tests by Type
% Cross-factor tests end with _Visual or _Saccade
% Cross-time tests contain _Vis_x_Sac or _Sac_x_Vis
test_names = {gen_tests_plan.test_name};
is_visual_epoch = endsWith(test_names, '_Visual');
is_saccade_epoch = endsWith(test_names, '_Saccade');
is_cross_epoch = contains(test_names, '_Vis_x_Sac') | ...
                 contains(test_names, '_Sac_x_Vis');

%% Define Factor Organization
% Columns = training factor, Rows = test factor
factors = {'Reward', 'Salience', 'Identity', 'Probability'};
n_factors = length(factors);

% Cross-epoch test types (for right panel)
cross_epoch_types = {'Vis_x_Sac', 'Sac_x_Vis'};

%% Identify Monkeys for Marker Styling
monkey_initials = cellfun(@(x) x(1), session_ids, 'UniformOutput', false);
[unique_monkeys, ~, ~] = unique(monkey_initials);
monkey_markers = {'o', 's', '^', 'd', 'p', 'h'};
if length(unique_monkeys) > length(monkey_markers)
    warning('More monkeys than defined markers; some will share markers.');
    marker_map = @(idx) monkey_markers{mod(idx-1, length(monkey_markers)) + 1};
else
    marker_map = @(idx) monkey_markers{idx};
end

%% Setup Figure - 90% of screen size
screen_size = get(0, 'ScreenSize');
fig_width = screen_size(3) * 0.9;
fig_height = screen_size(4) * 0.9;
fig_left = (screen_size(3) - fig_width) / 2;
fig_bottom = (screen_size(4) - fig_height) / 2;

fig = figure('Position', [fig_left, fig_bottom, fig_width, fig_height], ...
    'Color', 'w', 'MenuBar', 'None', 'ToolBar', 'None');

%% Define Viewport Regions (Pattern 3: Independent Column Groups)
% Layout: [Visual 4x4] [Saccade 4x4] [Cross-Epoch 4x2]
% Width allocation: ~38% Visual, ~38% Saccade, ~24% Cross-Epoch

% Common parameters
total_width = 0.88;
total_height = 0.82;
bottom_margin = 0.10;
left_margin_base = 0.08;
group_spacing = 0.04;  % Space between groups
title_space = 0.04;    % Space for group titles at top

% Calculate widths proportionally (4:4:2 ratio for columns)
visual_width = total_width * (4/10) - group_spacing/2;
saccade_width = total_width * (4/10) - group_spacing/2;
cross_width = total_width * (2/10);

% Calculate left margins for each group
visual_left = left_margin_base;
saccade_left = visual_left + visual_width + group_spacing;
cross_left = saccade_left + saccade_width + group_spacing;

% Viewport options for each group
opts_visual = {'Width', visual_width, 'Height', total_height - title_space, ...
    'LeftMargin', visual_left, 'BottomMargin', bottom_margin, ...
    'WidthSpacing', 0.008, 'HeightSpacing', 0.025};
opts_saccade = {'Width', saccade_width, 'Height', total_height - title_space, ...
    'LeftMargin', saccade_left, 'BottomMargin', bottom_margin, ...
    'WidthSpacing', 0.008, 'HeightSpacing', 0.025};
opts_cross = {'Width', cross_width, 'Height', total_height - title_space, ...
    'LeftMargin', cross_left, 'BottomMargin', bottom_margin, ...
    'WidthSpacing', 0.012, 'HeightSpacing', 0.025};

%% Create axes storage
h_axes_visual = gobjects(n_factors, n_factors);
h_axes_saccade = gobjects(n_factors, n_factors);
h_axes_cross = gobjects(n_factors, 2);

%% Plot Visual Epoch Group (4x4)
for i_row = 1:n_factors  % test factor (row)
    test_factor = factors{i_row};
    for i_col = 1:n_factors  % training factor (column)
        train_factor = factors{i_col};

        % Skip diagonal (self-test) and Salience×Identity
        if i_row == i_col
            continue;  % Diagonal: can't test on self
        end
        if (strcmp(train_factor, 'Salience') && strcmp(test_factor, 'Identity')) || ...
           (strcmp(train_factor, 'Identity') && strcmp(test_factor, 'Salience'))
            continue;  % Salience×Identity impossible
        end

        % Build test name: {TrainFactor}_x_{TestFactor}_Visual
        test_name = sprintf('%s_x_%s_Visual', train_factor, test_factor);

        % Find matching plan item
        plan_idx = find(strcmp(test_names, test_name));
        if isempty(plan_idx)
            continue;
        end
        plan_item = gen_tests_plan(plan_idx);

        % Create subplot (row-major indexing)
        subplot_idx = (i_row - 1) * n_factors + i_col;
        h_axes_visual(i_row, i_col) = mySubPlot([n_factors, n_factors, subplot_idx], ...
            opts_visual{:});

        % Plot data
        plot_single_panel(h_axes_visual(i_row, i_col), plan_item, ...
            aggregated_data, unique_monkeys, marker_map);

        % Add title (abbreviated)
        title(sprintf('%s→%s', train_factor(1), test_factor(1)), ...
            'FontSize', 8, 'FontWeight', 'normal');
    end
end

%% Plot Saccade Epoch Group (4x4)
for i_row = 1:n_factors  % test factor (row)
    test_factor = factors{i_row};
    for i_col = 1:n_factors  % training factor (column)
        train_factor = factors{i_col};

        % Skip diagonal and Salience×Identity
        if i_row == i_col
            continue;
        end
        if (strcmp(train_factor, 'Salience') && strcmp(test_factor, 'Identity')) || ...
           (strcmp(train_factor, 'Identity') && strcmp(test_factor, 'Salience'))
            continue;
        end

        % Build test name
        test_name = sprintf('%s_x_%s_Saccade', train_factor, test_factor);

        % Find matching plan item
        plan_idx = find(strcmp(test_names, test_name));
        if isempty(plan_idx)
            continue;
        end
        plan_item = gen_tests_plan(plan_idx);

        % Create subplot
        subplot_idx = (i_row - 1) * n_factors + i_col;
        h_axes_saccade(i_row, i_col) = mySubPlot([n_factors, n_factors, subplot_idx], ...
            opts_saccade{:});

        % Plot data
        plot_single_panel(h_axes_saccade(i_row, i_col), plan_item, ...
            aggregated_data, unique_monkeys, marker_map);

        % Add title
        title(sprintf('%s→%s', train_factor(1), test_factor(1)), ...
            'FontSize', 8, 'FontWeight', 'normal');
    end
end

%% Plot Cross-Epoch Group (4x2)
for i_row = 1:n_factors  % factor
    factor = factors{i_row};
    for i_col = 1:2  % Vis_x_Sac or Sac_x_Vis
        cross_type = cross_epoch_types{i_col};

        % Build test name
        test_name = sprintf('%s_%s', factor, cross_type);

        % Find matching plan item
        plan_idx = find(strcmp(test_names, test_name));
        if isempty(plan_idx)
            continue;
        end
        plan_item = gen_tests_plan(plan_idx);

        % Create subplot
        subplot_idx = (i_row - 1) * 2 + i_col;
        h_axes_cross(i_row, i_col) = mySubPlot([n_factors, 2, subplot_idx], ...
            opts_cross{:});

        % Plot data
        plot_single_panel(h_axes_cross(i_row, i_col), plan_item, ...
            aggregated_data, unique_monkeys, marker_map);

        % Add title
        if i_col == 1
            title_str = sprintf('%s V→S', factor(1));
        else
            title_str = sprintf('%s S→V', factor(1));
        end
        title(title_str, 'FontSize', 8, 'FontWeight', 'normal');
    end
end

%% Standardize All Axes
all_axes = [h_axes_visual(:); h_axes_saccade(:); h_axes_cross(:)];
valid_axes = all_axes(isgraphics(all_axes));

if isempty(valid_axes)
    warning('No valid axes created. Check data availability.');
    return;
end

set(valid_axes, 'TickDir', 'Out', 'Color', 'none', 'LineWidth', 0.5, ...
    'XColor', 'k', 'YColor', 'k', 'FontSize', 7, ...
    'XLim', [0.3 0.9], 'YLim', [0.3 0.9], ...
    'XTick', [0.4 0.6 0.8], 'YTick', [0.4 0.6 0.8]);

%% Remove Inner Tick Labels
% Visual group
for i_row = 1:n_factors
    for i_col = 1:n_factors
        if isgraphics(h_axes_visual(i_row, i_col))
            if i_row < n_factors  % Not bottom row
                set(h_axes_visual(i_row, i_col), 'XTickLabel', []);
            end
            if i_col > 1  % Not leftmost column
                set(h_axes_visual(i_row, i_col), 'YTickLabel', []);
            end
        end
    end
end

% Saccade group
for i_row = 1:n_factors
    for i_col = 1:n_factors
        if isgraphics(h_axes_saccade(i_row, i_col))
            if i_row < n_factors
                set(h_axes_saccade(i_row, i_col), 'XTickLabel', []);
            end
            if i_col > 1
                set(h_axes_saccade(i_row, i_col), 'YTickLabel', []);
            end
        end
    end
end

% Cross group
for i_row = 1:n_factors
    for i_col = 1:2
        if isgraphics(h_axes_cross(i_row, i_col))
            if i_row < n_factors
                set(h_axes_cross(i_row, i_col), 'XTickLabel', []);
            end
            if i_col > 1
                set(h_axes_cross(i_row, i_col), 'YTickLabel', []);
            end
        end
    end
end

%% Add Group Titles
title_y = bottom_margin + (total_height - title_space) + title_space/2;

% Visual Epoch title
annotation('textbox', [visual_left, title_y, visual_width, title_space], ...
    'String', 'Visual Epoch', 'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none');

% Saccade Epoch title
annotation('textbox', [saccade_left, title_y, saccade_width, title_space], ...
    'String', 'Saccade Epoch', 'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none');

% Cross-Epoch title
annotation('textbox', [cross_left, title_y, cross_width, title_space], ...
    'String', 'Cross-Epoch', 'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none');

%% Add Column Labels (Training Factor) for Visual and Saccade Groups
label_y = bottom_margin + (total_height - title_space) + 0.005;
for i_col = 1:n_factors
    % Visual group column labels
    col_center_vis = visual_left + (i_col - 0.5) * (visual_width / n_factors);
    annotation('textbox', [col_center_vis - 0.03, label_y, 0.06, 0.02], ...
        'String', factors{i_col}(1), 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'EdgeColor', 'none');

    % Saccade group column labels
    col_center_sac = saccade_left + (i_col - 0.5) * (saccade_width / n_factors);
    annotation('textbox', [col_center_sac - 0.03, label_y, 0.06, 0.02], ...
        'String', factors{i_col}(1), 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'EdgeColor', 'none');
end

%% Add Row Labels (Test Factor) for Visual Group
for i_row = 1:n_factors
    row_center = bottom_margin + (n_factors - i_row + 0.5) * ...
        ((total_height - title_space) / n_factors);
    annotation('textbox', [visual_left - 0.04, row_center - 0.015, 0.035, 0.03], ...
        'String', factors{i_row}(1), 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none');
end

%% Add Legend
legend_ax = axes('Position', [cross_left + cross_width - 0.06, ...
    bottom_margin + 0.02, 0.05, 0.08], 'Visible', 'off');
hold(legend_ax, 'on');
legend_handles = gobjects(1, length(unique_monkeys));
for i = 1:length(unique_monkeys)
    legend_handles(i) = plot(legend_ax, nan, nan, marker_map(i), ...
        'MarkerSize', 6, 'LineStyle', 'none', ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
end
legend(legend_ax, legend_handles, unique_monkeys, 'Location', 'best', ...
    'Box', 'off', 'FontSize', 8);

%% Add Global Axis Labels
% Y-axis label (left side)
annotation('textbox', [0.01, 0.3, 0.03, 0.4], ...
    'String', sprintf('Generalization Accuracy - %s', brain_area_name), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none');

% X-axis label (bottom)
annotation('textbox', [0.1, 0.01, 0.8, 0.04], ...
    'String', sprintf('Standard CV Accuracy - %s', brain_area_name), ...
    'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none');

%% Save Figure
figures_dir = fullfile(project_root, 'figures', 'decoding');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_decoding_generalization_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);
end

%% Helper Function: Plot Single Panel
function plot_single_panel(ax, plan_item, aggregated_data, unique_monkeys, marker_map)
    axes(ax);
    hold on;

    test_name = plan_item.test_name;

    % Check for generalization data
    if ~isfield(aggregated_data.population_decoding, test_name)
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', ...
            'FontSize', 7, 'Color', [0.5 0.5 0.5]);
        return;
    end
    gen_data = aggregated_data.population_decoding.(test_name);

    % Get standard CV test
    train_model_tag = plan_item.train_model_tag;
    standard_test_name = [train_model_tag, '_CV'];

    if ~isfield(aggregated_data.population_decoding, standard_test_name)
        text(0.5, 0.5, 'No CV', 'HorizontalAlignment', 'center', ...
            'FontSize', 7, 'Color', [0.5 0.5 0.5]);
        return;
    end
    standard_data = aggregated_data.population_decoding.(standard_test_name);

    % Create map for quick session lookup
    standard_map = containers.Map({standard_data.session_id}, ...
        1:length(standard_data));

    % Plot data for each session
    for i_session = 1:length(gen_data)
        session_id = gen_data(i_session).session_id;

        if ~isKey(standard_map, session_id) || isempty(session_id)
            continue;
        end

        std_idx = standard_map(session_id);

        x_val = standard_data(std_idx).accuracy;
        x_ci = standard_data(std_idx).accuracy_ci;
        y_val = gen_data(i_session).accuracy;
        y_ci = gen_data(i_session).accuracy_ci;

        % Determine marker by monkey
        monkey_initial = session_id(1);
        monkey_idx = find(strcmp(unique_monkeys, monkey_initial));
        if isempty(monkey_idx), continue; end
        marker_style = marker_map(monkey_idx);

        % Plot with error bars
        errorbar(x_val, y_val, y_val - y_ci(1), y_ci(2) - y_val, ...
            x_val - x_ci(1), x_ci(2) - x_val, ...
            marker_style, 'MarkerSize', 5, 'LineWidth', 0.5, ...
            'CapSize', 0, 'MarkerEdgeColor', 'none', ...
            'MarkerFaceColor', 'k', 'Color', [0.7 0.7 0.7]);
    end

    % Add reference lines
    plot([0.3 0.9], [0.5 0.5], 'k:', 'LineWidth', 0.5);
    plot([0.5 0.5], [0.3 0.9], 'k:', 'LineWidth', 0.5);
end