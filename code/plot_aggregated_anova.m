%% plot_aggregated_anova.m
%
% Generates a summary figure visualizing the proportion of significant
% neurons over time for each ANOVA model term across different alignment
% events. It compares SC and SNc populations.
%
% This function is designed to work with a nested data structure where
% results are organized first by alignment event, then by p-value type:
% aggregated_data.anova_results.(alignment_event).(p_value_field)
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-12
%
function plot_aggregated_anova(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
% Define the project root, ensure the figures directory exists, and add the
% 'utils' directory to the path so that helper functions can be found.
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Dynamically Discover Alignment Events and ANOVA Terms
% The function discovers the structure of the analysis from the data itself,
% making it robust to different analysis plans.

% Get alignment events (columns of the plot)
alignment_events = fieldnames(aggregated_sc_data.anova_results);
n_cols = numel(alignment_events);

% Get a unique superset of all p-value fields (rows of the plot)
try
all_p_fields = {};
for i_event = 1:n_cols
    event_name = alignment_events{i_event};
    all_p_fields = [all_p_fields; fieldnames(aggregated_sc_data.anova_results.(event_name))];
end
catch me
    keyboard
end

p_value_fields = unique(all_p_fields);
n_rows = numel(p_value_fields);

if n_rows == 0 || n_cols == 0
    warning('No ANOVA results found in the aggregated data. Cannot plot.');
    return;
end

%% Setup Figure
fig = figure('Position', [100, 100, 1200, 900], 'Color', 'w');
sc_color = [0, 0.4470, 0.7410];  % Blue
snc_color = [0.8500, 0.3250, 0.0980]; % Orange
h_axes = gobjects(n_rows, n_cols); % Use a 2D array for axes handles

%% Nested Plotting Loop (Row = ANOVA term, Col = Alignment Event)
for i_row = 1:n_rows
    p_value_name = p_value_fields{i_row};

    for i_col = 1:n_cols
        event_name = alignment_events{i_col};

        % --- Select Subplot ---
        plot_idx = (i_row - 1) * n_cols + i_col;
        h_axes(i_row, i_col) = mySubPlot([n_rows, n_cols, plot_idx]);
        hold on;

        % --- Data Processing and Plotting ---
        % Check if this specific p-value exists for this event
        if isfield(aggregated_sc_data.anova_results.(event_name), ...
                p_value_name)
            
            % Process and Plot SC Data
            p_values_sc = aggregated_sc_data.anova_results.( ...
                event_name).(p_value_name);
            prop_sig_sc = mean(p_values_sc < 0.05, 1, 'omitnan');
            n_time_bins = size(prop_sig_sc, 2);
            time_vector = 1:n_time_bins; % Generic time bins
            
            plot(time_vector, prop_sig_sc, 'Color', sc_color, ...
                'LineWidth', 2);

            % Process and Plot SNc Data
            p_values_snc = aggregated_snc_data.anova_results.( ...
                event_name).(p_value_name);
            prop_sig_snc = mean(p_values_snc < 0.05, 1, 'omitnan');

            plot(time_vector, prop_sig_snc, 'Color', snc_color, ...
                'LineWidth', 2);
        else
            % If data doesn't exist, display a note on the plot
            text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center');
            axis off;
        end
        
        box off;
    end
end

%% De-clutter Axes and Add Labels
for i_row = 1:n_rows
    for i_col = 1:n_cols
        ax = h_axes(i_row, i_col);
        
        % Add Column Titles (Alignment Events) to the top row
        if i_row == 1
            title(ax, strrep(alignment_events{i_col}, '_', ' '), ...
                'FontWeight', 'normal');
        end

        % Add Row Labels (ANOVA Terms) to the first column
        if i_col == 1
            clean_label = strrep(p_value_fields{i_row}, 'p_', '');
            clean_label = strrep(clean_label, '_', ' x ');
            ylabel(ax, clean_label, 'Interpreter', 'none');
        end

        % Remove X-tick labels from all but the bottom row
        if i_row < n_rows
           set(ax, 'XTickLabel', []);
        else
           xlabel(ax, 'Time Bins');
        end
        
        % Remove Y-tick labels from all but the first column
        if i_col > 1
            set(ax, 'YTickLabel', []);
        end
    end
end

%% Final Figure Touches
% Link all Y-axes to ensure consistent scaling for comparison
linkaxes(h_axes(:), 'y');
ylim(h_axes(1,1), [0, 0.5]); % Set a consistent y-limit for all plots

% Create a single, clear legend for the entire figure
lgd = legend(h_axes(1,1), {'SC', 'SNc'}, 'Location', 'northeast');
lgd.Box = 'off';

% Add an overarching title for the entire figure
sgtitle(['Proportion of Significant Neurons by ANOVA Term and ' ...
    '    Alignment Event'], 'FontWeight', 'bold');

% Save figure:
fig_filename = fullfile(figures_dir, 'aggregated_anova.pdf');
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end