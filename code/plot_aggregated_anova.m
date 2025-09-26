%% plot_aggregated_anova.m
%
%   Generates a summary figure visualizing time-resolved ANOVA results.
%
%   This function plots the proportion of neurons with a significant effect
%   (p < 0.05) for each term in the ANOVA models. The figure is a grid
%   with ANOVA terms as rows and alignment events as columns. Each subplot
%   shows the mean proportion across all sessions (bold trace) overlaid on
%   the traces from individual sessions (faint traces).
%
% INPUTS:
%
%   aggregated_data - Struct with aggregated data for one brain area,
%                     containing the `anova_results` field.
%
%   brain_area_name - String with the name of the brain area (e.g., 'SC').
%
%   analysis_plan   - The analysis plan struct, which is the source of
%                     truth for event names and model terms.
%
% Author: Jules
% Date: 2025-09-26 (Refactored for style and capitalization fix)
%
function plot_aggregated_anova(aggregated_data, brain_area_name, ...
    analysis_plan)

%% Setup
project_root = fullfile(findOneDrive, 'Code', ...
    '4factors-analysis-pipeline');
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
    
    % Derive a label from the model name (e.g., 'Image' or 'Bullseye')
    if contains(model_name, 'image')
        model_label = 'Image';
    elseif contains(model_name, 'bullseye')
        model_label = 'Bullseye';
    else
        model_label = 'Model';
    end

    % Dynamically generate all possible ANOVA terms up to pairwise
    % interactions from the plan's factors
    factors = anova_plan(i_model).factors;
    all_terms = {};
    for k = 1:2
        combs = nchoosek(factors, k);
        for i = 1:size(combs, 1)
            all_terms{end+1} = strjoin(combs(i,:), ':');
        end
    end

    for i_term = 1:length(all_terms)
        term = all_terms{i_term};
        
        % Capitalize term parts to match the field names in the data struct
        term_parts = strsplit(term, ':');
        capitalized_parts = cellfun(@(c) [upper(c(1)), c(2:end)], ...
            term_parts, 'UniformOutput', false);
        capitalized_term_for_field = strjoin(capitalized_parts, '_');
        p_field = ['p_', matlab.lang.makeValidName(...
            capitalized_term_for_field)];
        
        label = sprintf('%s (%s)', term, model_label);
        plot_layout = [plot_layout; {label, p_field, model_name}];
    end
    
    % Add a separator row between models if there are multiple
    if i_model < length(anova_plan)
        plot_layout = [plot_layout; {'---', '', ''}];
    end
end

n_rows = size(plot_layout, 1);
event_names = analysis_plan.events;
n_cols = length(event_names);

fig = figure('Position', [1, 1, 110 * n_cols, 100 * n_rows], ...
    'Color', 'w', 'MenuBar', 'None', 'ToolBar', 'None');
h_axes = gobjects(n_rows, n_cols);
colors = richColors;
colors = colors([6, 12], :);

%% Plotting Loop
for i_row = 1:n_rows
    label_name = plot_layout{i_row, 1};
    p_field = plot_layout{i_row, 2};
    analysis_name = plot_layout{i_row, 3};

    % Skip separator rows in the layout
    if strcmp(label_name, '---')
        continue;
    end

    for i_col = 1:n_cols
        event_name = event_names{i_col};
        
        ax_idx = (i_row - 1) * n_cols + i_col;
        ax = mySubPlot([n_rows, n_cols, ax_idx], 'Width', 0.875, ...
            'Height', 0.9, 'LeftMargin', 0.1);
        h_axes(i_row, i_col) = ax;
        hold on;
        
        % Check if the required data exists
        if ~isfield(aggregated_data.anova_results, analysis_name) || ...
           ~isfield(aggregated_data.anova_results.(analysis_name), ...
           event_name)
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off; continue;
        end

        session_results = ...
            aggregated_data.anova_results.(analysis_name).(event_name);
        n_sessions = length(session_results);
        all_proportions = [];
        time_vector = [];

        % Aggregate proportions from all sessions for this subplot
        for i_session = 1:n_sessions
            if isfield(session_results(i_session), p_field) && ...
                    ~isempty(session_results(i_session).(p_field))

                p_values = session_results(i_session).(p_field);
                
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

        % Plot individual session traces with transparency
        for i_session = 1:size(all_proportions, 1)
            plot(ax, time_vector, all_proportions(i_session, :), ...
                'Color', [colors(1,:), 0.2]);
        end
        
        % Plot the mean trace across sessions
        mean_trace = mean(all_proportions, 1, 'omitnan');
        plot(ax, time_vector, mean_trace, 'Color', colors(1,:), ...
            'LineWidth', 2);

        xlim(ax, [time_vector(1), time_vector(end)]);
        box off;
    end
end

%% Final Figure Formatting
% Use outerLims to set common Y-limits for each row of plots
for i_row = 1:n_rows
    % Skip separator rows
    if strcmp(plot_layout{i_row, 1}, '---')
        continue;
    end
    
    % Get all valid axes handles for the current row
    row_axes = h_axes(i_row, isgraphics(h_axes(i_row, :)));
    if isempty(row_axes)
        continue;
    end
    
    % Get the outer limits for this row and set a common y-axis
    [~, yLims] = outerLims(row_axes);
    set(row_axes, 'YLim', [0, max(0.1, yLims(2) * 1.1)]);
end

% Link x-axes for all valid plots
valid_axes_mask = ~strcmp(plot_layout(:,1),'---');
linkaxes(h_axes(valid_axes_mask, :), 'x');

% De-clutter axes by removing redundant labels
for i_row = 1:n_rows
    if strcmp(plot_layout{i_row, 1}, '---')
        continue;
    end
    for i_col = 1:n_cols
        ax = h_axes(i_row, i_col);
        if ~isgraphics(ax)
            continue;
        end

        % Remove Y-tick labels for all but the first column of plots
        if i_col > 1
            set(ax, 'YTickLabel', []);
        end

        % Set X-label only for the bottom-most plots in the grid
        is_last_row = (i_row == n_rows);
        is_last_visible_row = is_last_row || ...
            strcmp(plot_layout{i_row + 1, 1}, '---');
        % Column Labeling (Bottom-most row)
        if is_last_visible_row
             xlabel_str = sprintf('Time from %s (s)', ...
                 strrep(event_names{i_col}, '_', ' '));
             xlabel(ax, xlabel_str);
        else
            set(ax, 'XTickLabel', []);
            xlabel(ax, '');
        end
    end
end

%% Add Overarching Y-Labels for Model Blocks
% --- Add Overarching Y-Labels for Model Blocks ---

% Find the row indices corresponding to each ANOVA model
model_names = plot_layout(:, 3);
image_rows = find(strcmp(model_names, 'anova_imagetrials'));
bullseye_rows = find(strcmp(model_names, 'anova_bullseyetrials'));

% --- Add Label for Image Trials Block ---
if ~isempty(image_rows)
    % Get handles for the top-left and bottom-left axes in this block
    top_left_ax = h_axes(image_rows(1), 1);
    bottom_left_ax = h_axes(image_rows(end), 1);

    % Get their positions to calculate the bounding box
    top_pos = get(top_left_ax, 'Position');
    bottom_pos = get(bottom_left_ax, 'Position');

    % Calculate position for the new invisible axes
    label_ax_pos = bottom_pos;
    label_ax_pos(4) = (top_pos(2) + top_pos(4)) - bottom_pos(2);

    % Create the invisible axes and apply the ylabel
    label_ax = axes('Position', label_ax_pos, 'Visible', 'off');
    ylabel_str = sprintf('Image Trials - %s', brain_area_name);
    ylabel(label_ax, ylabel_str, 'Visible', 'on', ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
end

% --- Add Label for Bullseye Trials Block ---
if ~isempty(bullseye_rows)
    % Get handles for the top-left and bottom-left axes in this block
    top_left_ax = h_axes(bullseye_rows(1), 1);
    bottom_left_ax = h_axes(bullseye_rows(end), 1);

    % Get their positions to calculate the bounding box
    top_pos = get(top_left_ax, 'Position');
    bottom_pos = get(bottom_left_ax, 'Position');

    % Calculate position for the new invisible axes
    label_ax_pos = bottom_pos;
    label_ax_pos(4) = (top_pos(2) + top_pos(4)) - bottom_pos(2);

    % Create the invisible axes and apply the ylabel
    label_ax = axes('Position', label_ax_pos, 'Visible', 'off');
    ylabel_str = sprintf('Bullseye Trials - %s', brain_area_name);
    ylabel(label_ax, ylabel_str, 'Visible', 'on', ...
        'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
end

%% Save Figure
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_anova_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end