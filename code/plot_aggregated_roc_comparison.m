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

%% Dynamically Discover Grid Layout
event_names = fieldnames(aggregated_data.roc_comparison);
n_events = length(event_names);

% Collect all unique comparison names from all events
all_comparison_names = {};
for i_event = 1:n_events
    event_name = event_names{i_event};
    all_comparison_names = [all_comparison_names; fieldnames(aggregated_data.roc_comparison.(event_name))];
end
comparison_names = unique(all_comparison_names, 'stable');
n_comparisons = length(comparison_names);

if n_comparisons == 0 || n_events == 0
    warning('plot_aggregated_roc_comparison:no_data_to_plot', ...
        'No valid comparisons or events found to plot.');
    return;
end

%% Figure and Plotting Setup
fig = figure('Position', [100, 100, 300 * n_events, 250 * n_comparisons], 'Color', 'w');
h_axes = gobjects(n_comparisons, n_events);

colors = richColors;
posColor = colors(7,:); % Color for cond2 > cond1
negColor = colors(10,:); % Color for cond1 > cond2

%% Main Plotting Loop (Rows = Comparisons, Columns = Events)
for i_row = 1:n_comparisons
    comp_name = comparison_names{i_row};

    for i_col = 1:n_events
        event_name = event_names{i_col};

        % --- Subplot Creation ---
        h_axes(i_row, i_col) = mySubPlot([n_comparisons, n_events, i_col + (i_row-1)*n_events]);

        % --- Data Retrieval ---
        % Check if this specific comparison exists for this event
        if ~isfield(aggregated_data.roc_comparison.(event_name), comp_name)
            title(sprintf('%s (N/A)', strrep(comp_name, '_', ' ')));
            set(gca, 'Color', [0.95 0.95 0.95]); % Indicate missing data
            continue;
        end

        session_data = aggregated_data.roc_comparison.(event_name).(comp_name);

        % --- Data Aggregation from Struct Array ---
        if isempty(session_data)
            title(sprintf('%s (No Data)', strrep(comp_name, '_', ' ')));
            continue;
        end

        all_sig_matrices = arrayfun(@(s) s.sig, session_data, 'UniformOutput', false);
        sig_matrix = cat(1, all_sig_matrices{:});
        time_vector = session_data(1).time_vector;

        % --- Data Processing ---
        if isempty(sig_matrix)
            prop_cond2_pref = zeros(1, length(time_vector));
            prop_cond1_pref = zeros(1, length(time_vector));
        else
            prop_cond2_pref = mean(sig_matrix == 1, 1, 'omitnan');
            prop_cond1_pref = -mean(sig_matrix == -1, 1, 'omitnan');
        end

        % --- Plotting Data ---
        hold on;
        h1 = barStairsFill(time_vector, zeros(size(prop_cond2_pref)), prop_cond2_pref);
        delete(h1(2)); % Delete the line from barStairsFill
        set(h1(1), 'FaceColor', posColor, 'EdgeColor', 'none');
        set(h1(3), 'Color', posColor);

        h2 = barStairsFill(time_vector, zeros(size(prop_cond1_pref)), prop_cond1_pref);
        delete(h2(2)); % Delete the line from barStairsFill
        set(h2(1), 'FaceColor', negColor, 'EdgeColor', 'none');
        set(h2(3), 'Color', negColor);

        xlim([time_vector(1), time_vector(end)]);
        line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
        line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
    end
end

%% Final Figure Formatting
% --- Labeling and Axis Standardization ---
for i_row = 1:n_comparisons
    % Standardize Y-axis limits for the current row
    valid_axes_in_row = h_axes(i_row, isgraphics(h_axes(i_row, :)));
    if ~isempty(valid_axes_in_row)
        [~, yLims] = outerLims(valid_axes_in_row);
        max_abs_y = max(abs(yLims));
        if max_abs_y == 0; max_abs_y = 1; end % Avoid flat limits
        set(valid_axes_in_row, 'YLim', [-max_abs_y, max_abs_y]);
    end

    for i_col = 1:n_events
        h_ax = h_axes(i_row, i_col);
        if ~isgraphics(h_ax); continue; end % Skip empty/invalid subplots

        % --- X-axis labels (bottom row only) ---
        if i_row == n_comparisons
            xlabel_str = sprintf('Time from %s (s)', strrep(event_names{i_col}, '_', ' '));
            xlabel(h_ax, xlabel_str);
        else
            set(h_ax, 'XTickLabel', []);
        end

        % --- Y-axis labels (first column only) ---
        if i_col == 1
            comp_name_str = strrep(comparison_names{i_row}, '_', ' ');
            ylabel(h_ax, comp_name_str);
        else
            set(h_ax, 'YTickLabel', []);
        end
    end
end

% --- General Formatting for all valid axes ---
all_valid_axes = h_axes(isgraphics(h_axes));
set(all_valid_axes, 'TickDir', 'Out', 'Color', 'none', ...
    'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

%% Save Figure
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_roc_comparison_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
