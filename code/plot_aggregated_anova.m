%% plot_aggregated_anova.m
%
%   Generates a summary figure visualizing time-resolved ANOVA results,
%   adhering to the specifications in `docs/plotting_requirements.md`.
%
%   It plots the proportion of significant neurons (p < 0.05) for each
%   session individually, with the cross-session average trace overlaid.
%   The figure is structured with each ANOVA effect term as a row and
%   each alignment event as a column.
%
% INPUTS:
%
%   aggregated_data - A struct containing aggregated data for one brain area.
%                     This script requires the `anova_results` field, structured as:
%                     .(analysis_name).(event_name) = [1 x nSessions] struct array
%
%                     Each struct in the array must contain fields for the p-value
%                     of each model term (e.g., .p_reward) and a .time_vector.
%
%   brain_area_name - A string with the name of the brain area (e.g., 'SC').
%   analysis_plan   - The analysis plan struct, typically loaded from the same
%                     .mat file as aggregated_data. This is the source of truth
%                     for event names and model terms.
%
% Author: Jules
% Date: 2025-09-24
%
function plot_aggregated_anova(aggregated_data, brain_area_name, analysis_plan)

%% Setup
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

if ~isfield(aggregated_data, 'anova_results')
    warning('plot_aggregated_anova:no_data', 'No ANOVA data found.');
    return;
end

%% Dynamically Determine Plot Layout from Analysis Plan
anova_plan = analysis_plan.anova_plan;
plot_layout = {};

for i_model = 1:length(anova_plan)
    model_name = anova_plan(i_model).name;
    
    % --- Start of Corrected Logic ---
    
    % 1. Derive a label from the model name
    if contains(model_name, 'image')
        model_label = 'Image';
    elseif contains(model_name, 'bullseye')
        model_label = 'Bullseye';
    else
        model_label = 'Model';
    end

    % 2. Dynamically generate all possible ANOVA terms from the .factors field
    factors = anova_plan(i_model).factors;
    all_terms = {};
    for k = 1:length(factors)
        combs = nchoosek(factors, k);
        for i = 1:size(combs, 1)
            all_terms{end+1} = strjoin(combs(i,:), ':');
        end
    end
    
    % --- End of Corrected Logic ---

    for i_term = 1:length(all_terms)
        term = all_terms{i_term};
        % Create the p-value field name, replacing ':' with '_' for struct 
        % compatibility
        p_field = ['p_', matlab.lang.makeValidName(strrep(term, ':', '_'))];
        label = sprintf('%s (%s)', term, model_label);
        plot_layout = [plot_layout; {label, p_field, model_name}];
    end
    
    % Add a separator between models if there are multiple
    if i_model < length(anova_plan)
        plot_layout = [plot_layout; {'---', '', ''}];
    end
end

n_rows = size(plot_layout, 1);
event_names = analysis_plan.events;
n_cols = length(event_names);

fig = figure('Position', [100, 100, 300 * n_cols, 150 * n_rows], ...
    'Color', 'w');
h_axes = gobjects(n_rows, n_cols);
colors = richColors;
colors = colors([6, 12], :);

%% Plotting Loop
for i_row = 1:n_rows
    label_name = plot_layout{i_row, 1};
    p_field = plot_layout{i_row, 2};
    analysis_name = plot_layout{i_row, 3};

    if strcmp(label_name, '---')
        continue; % Skip separator row
    end

    for i_col = 1:n_cols
        event_name = event_names{i_col};
        ax = mySubPlot([n_rows, n_cols, (i_row - 1) * n_cols + i_col]);
        h_axes(i_row, i_col) = ax;
        hold on;

        if ~isfield(aggregated_data.anova_results, analysis_name) || ...
           ~isfield(aggregated_data.anova_results.(analysis_name), ...
           event_name)
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off; continue;
        end

        session_results = aggregated_data.anova_results.( ...
            analysis_name).(event_name);
        n_sessions = length(session_results);
        all_proportions = [];
        time_vector = [];

        % Aggregate proportions from all sessions
        for i_session = 1:n_sessions
            if isfield(session_results(i_session), p_field)

                p_values = session_results(i_session).(p_field);
                if isempty(p_values), continue; end

                if isempty(time_vector)
                    time_vector = session_results(i_session).time_vector;
                end

                session_proportion = mean(p_values < 0.05, 1, 'omitnan');
                all_proportions = [all_proportions; session_proportion];
            end
        end

        if isempty(all_proportions)
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off; continue;
        end

        % Plot individual session traces
        for i_session = 1:size(all_proportions, 1)
            plot(ax, time_vector, all_proportions(i_session, :), 'Color', [colors(1,:), 0.2]);
        end
        
        % Plot average trace
        mean_trace = mean(all_proportions, 1, 'omitnan');
        plot(ax, time_vector, mean_trace, 'Color', colors(1,:), 'LineWidth', 2);

        xlim(ax, [time_vector(1), time_vector(end)]);
        ylim(ax, [0, max(1, max(mean_trace) * 1.2)]); % Ensure Y-axis starts at 0
        box off;
    end
end

%% Final Figure Formatting
% Find a good y-limit across all related plots
y_max = 0;
for i_ax = 1:numel(h_axes)
    if isgraphics(h_axes(i_ax))
        current_ylim = get(h_axes(i_ax), 'YLim');
        if current_ylim(2) > y_max
            y_max = current_ylim(2);
        end
    end
end
if y_max == 0, y_max = 1; end % Default if no data

valid_axes_mask = ~strcmp(plot_layout(:,1),'---');
linkaxes(h_axes(valid_axes_mask, :), 'x');
set(h_axes(isgraphics(h_axes)), 'YLim', [0, y_max]);

for i_row = 1:n_rows
    if strcmp(plot_layout{i_row, 1}, '---'), continue; end
    for i_col = 1:n_cols
        ax = h_axes(i_row, i_col);
        if ~isgraphics(ax), continue; end

        % Set title only for the top row of plots
        if i_row == 1
            title(ax, strrep(event_names{i_col}, '_', ' '));
        end

        % Set Y-label only for the first column
        if i_col == 1
            ylabel(ax, plot_layout{i_row, 1}, 'Interpreter', 'none');
        else
            set(ax, 'YTickLabel', []);
        end

        % Set X-label only for the bottom-most plots in the grid
        is_last_visible_row = (i_row == n_rows) || strcmp(plot_layout{i_row+1, 1}, '---');
        if is_last_visible_row
            xlabel(ax, 'Time (s)');
        else
            set(ax, 'XTickLabel', []);
        end
    end
end

sgtitle(sprintf('Proportion of Significant Neurons (p < 0.05) for %s', brain_area_name), 'FontWeight', 'bold', 'Interpreter', 'none');

%% Save Figure
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_anova_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end