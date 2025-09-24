%% plot_aggregated_decoding.m
%
%   Generates a comprehensive summary figure visualizing population decoding
%   results. The output adheres to the specifications in
%   `docs/plotting_requirements.md`.
%
%   The figure is a grid of scatter plots, where each subplot represents a
%   specific decoding test defined in the analysis plan. Each point within a
%   plot corresponds to a single session.
%
% INPUTS:
%
%   aggregated_data - A struct containing aggregated data for one brain area.
%                     This script requires the `population_decoding` field,
%                     which is structured as:
%
%                     .population_decoding.(test_name) = [1 x nSessions] struct array
%
%                     Each struct in the array must contain:
%                       - .accuracy:    Mean classification accuracy.
%                       - .accuracy_ci: 95% confidence interval of the accuracy.
%                       - .session_id:  The unique identifier for the session.
%
%   brain_area_name - A string with the name of the brain area (e.g., 'SC').
%   session_ids     - A struct containing session IDs for each brain area.
%                     Used to identify monkeys for marker styling.
%
% Author: Jules
% Date: 2025-09-24
%

function plot_aggregated_decoding(aggregated_data, brain_area_name, session_ids)
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

%% Dynamically Discover Tests
all_tests = fieldnames(aggregated_data.population_decoding);
% Filter out the standard cross-validated tests to get the generalization tests
gen_tests = all_tests(~endsWith(all_tests, '_CV'));
n_tests = length(gen_tests);

if n_tests == 0
    disp('No generalization decoding tests found to plot.');
    return;
end

%% Setup Figure and Subplot Layout
n_cols = ceil(sqrt(n_tests));
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
    test_name = gen_tests{i_test};
    h_axes(i_test) = mySubPlot([n_rows, n_cols, i_test]);
    hold on;

    % --- Identify Standard and Generalization Data ---
    % The generalization data is for the current test
    gen_data = aggregated_data.population_decoding.(test_name);

    % The standard (x-axis) data comes from the corresponding CV test
    parts = strsplit(test_name, '_vs_');
    standard_test_name = [parts{1}, '_CV'];
    if ~isfield(aggregated_data.population_decoding, standard_test_name)
        warning('plot_aggregated_decoding:missing_cv', ...
            'Could not find standard CV test "%s" for generalization test "%s". Skipping.', ...
            standard_test_name, test_name);
        title(h_axes(i_test), {strrep(test_name, '_', ' '), '(Missing CV data)'});
        continue;
    end
    standard_data = aggregated_data.population_decoding.(standard_test_name);

    % Create maps for quick lookup of session data
    standard_map = containers.Map({standard_data.session_id}, 1:length(standard_data));

    % --- Plot Data for Each Session ---
    for i_session = 1:length(gen_data)
        session_id = gen_data(i_session).session_id;

        if ~isKey(standard_map, session_id)
            continue; % Skip if this session doesn't have a standard counterpart
        end

        std_idx = standard_map(session_id);

        % X-axis: Standard accuracy
        x_val = standard_data(std_idx).accuracy;
        x_ci  = standard_data(std_idx).accuracy_ci;

        % Y-axis: Generalization accuracy
        y_val = gen_data(i_session).accuracy;
        y_ci  = gen_data(i_session).accuracy_ci;

        % Determine marker for the monkey
        monkey_initial = session_id(1);
        monkey_idx = find(strcmp(unique_monkeys, monkey_initial));
        marker_style = marker_map(monkey_idx);

        % Plot point with X and Y error bars
        errorbar(x_val, y_val, ...
            y_val - y_ci(1), y_ci(2) - y_val, ... % y error
            x_val - x_ci(1), x_ci(2) - x_val, ... % x error
            marker_style, 'MarkerSize', 6, 'LineWidth', 1, 'CapSize', 0, 'Color', 'k');
    end

    % --- Formatting for each subplot ---
    axis([0 1 0 1]);
    plot([0 1], [0 1], 'k:', 'LineWidth', 0.5); % Diagonal line
    axis square;
    title(strrep(test_name, '_', ' '), 'Interpreter', 'none');
end

%% Figure Cleanup and Final Touches
set(h_axes, 'TickDir', 'Out', 'Color', 'none', 'LineWidth', 1, ...
    'XColor', 'k', 'YColor', 'k');

% Add a legend for the monkey markers
legend_handles = gobjects(1, length(unique_monkeys));
for i = 1:length(unique_monkeys)
    legend_handles(i) = plot(nan, nan, marker_map(i), ...
        'MarkerSize', 8, 'LineStyle', 'none', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
end
legend(legend_handles, unique_monkeys, 'Location', 'bestoutside', 'Box', 'off');

% Add shared labels
big_ax = axes('Position', [0.05 0.05 0.8 0.9], 'Visible', 'off');
ylabel(big_ax, 'Generalization Accuracy', 'Visible', 'on', 'FontSize', 12, 'FontWeight', 'bold');
xlabel(big_ax, 'Standard Accuracy', 'Visible', 'on', 'FontSize', 12, 'FontWeight', 'bold');

% Add Title and Save
sgtitle_str = sprintf('Decoding Generalization in %s', brain_area_name);
sgtitle(sgtitle_str, 'FontSize', 16, 'FontWeight', 'bold');

% Save the figure
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_decoding_generalization_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
