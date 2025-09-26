%% plot_aggregated_baseline_comparison.m
%
%   Generates a publication-quality summary figure for a single brain area,
%   showing the proportion of neurons with significant firing rate
%   increases or decreases compared to a baseline period. The output adheres
%   to the specifications in `docs/plotting_requirements.md`.
%
%   A single consolidated figure is generated with:
%       - Rows: Conditions (e.g., High Reward, Low Reward)
%       - Columns: Alignment Events (e.g., fixOn, targetOn)
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
% Date: 2025-09-26 (Refactored)
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

%% Dynamically Discover and Organize Plotting Grid
event_names = fieldnames(aggregated_data.baseline_comparison);
n_events = length(event_names);

% Discover all unique condition names across all events
all_condition_names = {};
for i_event = 1:n_events
    event_name = event_names{i_event};
    all_condition_names = [all_condition_names; fieldnames(aggregated_data.baseline_comparison.(event_name))];
end
condition_names = unique(all_condition_names, 'stable');
n_conditions = length(condition_names);

if n_conditions == 0 || n_events == 0
    warning('plot_aggregated_baseline_comparison:no_conditions_or_events', 'No conditions or events to plot.');
    return;
end

%% Pre-calculate Y-limits for each row to ensure consistent scaling
row_y_max = zeros(n_conditions, 1);
for i_cond = 1:n_conditions
    cond_name = condition_names{i_cond};
    max_val_for_row = 0;
    for i_event = 1:n_events
        event_name = event_names{i_event};
        if isfield(aggregated_data.baseline_comparison.(event_name), cond_name)
            session_data = aggregated_data.baseline_comparison.(event_name).(cond_name);
            if ~isempty(session_data)
                all_sig_matrices = arrayfun(@(s) s.sig, session_data, 'UniformOutput', false);
                sig_matrix = cat(1, all_sig_matrices{:});
                if ~isempty(sig_matrix)
                    prop_increase = mean(sig_matrix == 1, 1, 'omitnan');
                    prop_decrease = abs(mean(sig_matrix == -1, 1, 'omitnan'));
                    current_max = max([prop_increase, prop_decrease]);
                    if current_max > max_val_for_row
                        max_val_for_row = current_max;
                    end
                end
            end
        end
    end
    % Set a minimum limit and round up for clean axes
    row_y_max(i_cond) = ceil(max(max_val_for_row, 0.1) * 10) / 10;
end

%% Figure and Plotting Setup
fig_width = 250 * n_events;
fig_height = 200 * n_conditions;
fig = figure('Position', [100, 100, fig_width, fig_height], 'Color', 'w');

h_axes = gobjects(n_conditions, n_events);

colors = richColors();
plot_color = colors(1,:);
pos_alpha = 0.6; % Face alpha for increase
neg_alpha = 0.3; % Face alpha for decrease

%% Main Plotting Loop (Rows = Conditions, Columns = Events)
for i_cond = 1:n_conditions
    cond_name = condition_names{i_cond};
    y_lim_row = [-row_y_max(i_cond), row_y_max(i_cond)];

    for i_event = 1:n_events
        event_name = event_names{i_event};

        % Calculate subplot index (row-major order)
        plot_idx = (i_cond - 1) * n_events + i_event;
        h_axes(i_cond, i_event) = mySubPlot([n_conditions, n_events, plot_idx]);
        hold on;

        % Check if this condition exists for this event
        if ~isfield(aggregated_data.baseline_comparison.(event_name), cond_name)
            text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', 'FontSize', 9);
            axis off;
            continue;
        end

        session_data = aggregated_data.baseline_comparison.(event_name).(cond_name);

        % --- Data Aggregation ---
        if isempty(session_data)
            text(0.5, 0.5, 'No Data', 'HorizontalAlignment', 'center', 'FontSize', 9);
            axis off;
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
        h_inc = barStairsFill(time_vector, zeros(size(prop_increase)), prop_increase);
        set(h_inc(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', pos_alpha);
        delete(h_inc(2:3));

        h_dec = barStairsFill(time_vector, zeros(size(prop_decrease)), prop_decrease);
        set(h_dec(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', neg_alpha);
        delete(h_dec(2:3));

        % --- Formatting ---
        xlim([time_vector(1), time_vector(end)]);
        ylim(y_lim_row);
        line(xlim, [0, 0], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
        line([0, 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', ':');

    end % End of event loop (columns)
end % End of condition loop (rows)


%% Figure Cleanup and Final Touches
for i_cond = 1:n_conditions
    for i_event = 1:n_events
        ax = h_axes(i_cond, i_event);
        if ~isgraphics(ax, 'axes'), continue; end % Skip if 'No Data' plot

        set(ax, 'TickDir', 'out', 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

        % Row Labeling (Left-most column)
        if i_event == 1
            ylabel_str = strrep(condition_names{i_cond}, '_', ' ');
            ylabel(ax, ylabel_str, 'FontWeight', 'bold');
        else
            set(ax, 'YTickLabel', []);
            ylabel(ax, '');
        end

        % Column Labeling (Bottom-most row)
        if i_cond == n_conditions
             xlabel_str = sprintf('Time from %s (s)', strrep(event_names{i_event}, '_', ' '));
             xlabel(ax, xlabel_str);
        else
            set(ax, 'XTickLabel', []);
            xlabel(ax, '');
        end
    end
end

% Create a single, overarching Y-label for the entire figure
han = axes(fig, 'visible', 'off');
han.YLabel.Visible = 'on';
ylabel(han, 'Proportion of Neurons', 'FontWeight', 'bold');

% Add a main title
sgtitle(sprintf('Baseline Comparison for %s', strrep(brain_area_name, '_', ' ')), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

%% Save Figure
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_baseline_comparison_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end