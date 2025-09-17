%% plot_aggregated_anova.m
%
% Generates a summary figure for a single brain area, visualizing the
% proportion of significant neurons over time for each ANOVA model term
% across different alignment events.
%
% This function is designed to work with a nested data structure where
% results are organized first by alignment event, then by p-value type:
% aggregated_data.anova_results.(alignment_event).(p_value_field)
%
% INPUTS:
%   aggregated_data   - A struct containing aggregated data for the brain area.
%   brain_area_name   - A string with the name of the brain area (e.g., 'SC').
%
% Author: Jules
% Date: 2025-09-17
%
function plot_aggregated_anova(aggregated_data, brain_area_name)

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Input Validation
if ~isfield(aggregated_data, 'anova_results')
    warning('plot_aggregated_anova:no_data', ...
        'ANOVA results data is missing from the input struct.');
    return;
end

%% Dynamically Discover Alignment Events and ANOVA Terms
alignment_events = fieldnames(aggregated_data.anova_results);
n_cols = numel(alignment_events);

all_p_fields = {};
for i_event = 1:n_cols
    event_name = alignment_events{i_event};
    all_p_fields = [all_p_fields; fieldnames(aggregated_data.anova_results.(event_name))];
end
p_value_fields = unique(all_p_fields);
n_rows = numel(p_value_fields);

if n_rows == 0 || n_cols == 0
    warning('No ANOVA results found in the aggregated data. Cannot plot.');
    return;
end

%% Setup Figure
fig = figure('Position', [100, 100, 1200, 900], 'Color', 'w');
plot_color = [0, 0.4470, 0.7410];  % Blue
h_axes = gobjects(n_rows, n_cols);

%% Nested Plotting Loop (Row = ANOVA term, Col = Alignment Event)
for i_row = 1:n_rows
    p_value_name = p_value_fields{i_row};

    for i_col = 1:n_cols
        event_name = alignment_events{i_col};
        plot_idx = (i_row - 1) * n_cols + i_col;
        h_axes(i_row, i_col) = mySubPlot([n_rows, n_cols, plot_idx]);
        hold on;

        event_data = aggregated_data.anova_results.(event_name);

        if isfield(event_data, p_value_name)
            if isfield(event_data, 'time_vector')
                time_vector = event_data.time_vector;
            else
                n_time_bins = size(event_data.(p_value_name).p_values, 2);
                time_vector = 1:n_time_bins;
            end

            p_values = event_data.(p_value_name).p_values;
            prop_sig = mean(p_values < 0.05, 1, 'omitnan');
            
            plot(time_vector, prop_sig, 'Color', plot_color, 'LineWidth', 2);
            xlim(h_axes(i_row, i_col), [time_vector(1), time_vector(end)]);
        else
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
        
        if i_row == 1
            title(ax, strrep(alignment_events{i_col}, '_', ' '), 'FontWeight', 'normal');
        end

        if i_col == 1
            clean_label = strrep(p_value_fields{i_row}, 'p_', '');
            clean_label = strrep(clean_label, '_', ' x ');
            ylabel(ax, clean_label, 'Interpreter', 'none');
        end

        if i_row < n_rows
           set(ax, 'XTickLabel', []);
        else
           xlabel(ax, 'Time Bins');
        end
        
        if i_col > 1
            set(ax, 'YTickLabel', []);
        end
    end
end

%% Final Figure Touches
linkaxes(h_axes(:), 'y');
ylim(h_axes(1,1), [0, 0.5]);

sgtitle_str = sprintf('Proportion of Significant Neurons for %s by ANOVA Term', brain_area_name);
sgtitle(sgtitle_str, 'FontWeight', 'bold', 'Interpreter', 'none');

% Save figure:
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_anova_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end