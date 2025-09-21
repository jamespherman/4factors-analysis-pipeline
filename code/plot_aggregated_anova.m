%% plot_aggregated_anova.m
%
% Generates a summary figure visualizing time-resolved ANOVA results.
% It plots the proportion of significant neurons (p < 0.05) for each
% session individually, with the cross-session average trace overlaid.
% The figure is structured with each ANOVA effect term as a row and
% each alignment event as a column.
%
% INPUTS:
%   aggregated_data   - A struct with aggregated ANOVA results.
%   brain_area_name   - A string with the name of the brain area.
%
% Author: Jules
% Date: 2025-09-21
%
function plot_aggregated_anova(aggregated_data, brain_area_name)

%% Setup
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
addpath(fullfile(project_root, 'code', 'utils'));

if ~isfield(aggregated_data, 'anova_results')
    warning('plot_aggregated_anova:no_data', 'No ANOVA data found.');
    return;
end

% Define the 14 effect terms for the rows of the plot
effect_terms = { ...
    'Reward (Image)', 'p_reward', 'anova_imagetrials'; ...
    'Probability (Image)', 'p_probability', 'anova_imagetrials'; ...
    'Identity (Image)', 'p_identity', 'anova_imagetrials'; ...
    'R x P (Image)', 'p_reward_probability', 'anova_imagetrials'; ...
    'R x I (Image)', 'p_reward_identity', 'anova_imagetrials'; ...
    'P x I (Image)', 'p_probability_identity', 'anova_imagetrials'; ...
    'R x P x I (Image)', 'p_reward_probability_identity', 'anova_imagetrials'; ...
    '---', '', ''; ... % Separator
    'Reward (Bullseye)', 'p_reward', 'anova_bullseyetrials'; ...
    'Probability (Bullseye)', 'p_probability', 'anova_bullseyetrials'; ...
    'Saliency (Bullseye)', 'p_saliency', 'anova_bullseyetrials'; ...
    'R x P (Bullseye)', 'p_reward_probability', 'anova_bullseyetrials'; ...
    'R x S (Bullseye)', 'p_reward_saliency', 'anova_bullseyetrials'; ...
    'P x S (Bullseye)', 'p_probability_saliency', 'anova_bullseyetrials'; ...
    'R x P x S (Bullseye)', 'p_reward_probability_saliency', 'anova_bullseyetrials' ...
};
n_rows = size(effect_terms, 1);
[~, analysis_plan] = define_task_conditions();
n_cols = numel(analysis_plan.alignment_events); % nEvents

fig = figure('Position', [100, 100, 300 * n_cols, 100 * n_rows], 'Color', 'w');
h_axes = gobjects(n_rows, n_cols);
colors = richColors(2);

%% Plotting Loop
for i_row = 1:n_rows
    label_name = effect_terms{i_row, 1};
    p_field = effect_terms{i_row, 2};
    analysis_name = effect_terms{i_row, 3};

    if strcmp(label_name, '---')
        continue; % Skip separator row
    end

    for i_col = 1:n_cols
        event_name = analysis_plan.alignment_events{i_col};
        ax = mySubPlot([n_rows, n_cols, (i_row - 1) * n_cols + i_col]);
        h_axes(i_row, i_col) = ax;
        hold on;

        session_results = aggregated_data.anova_results.(analysis_name);
        n_sessions = length(session_results);
        all_proportions = [];
        time_vector = [];

        % Aggregate proportions from all sessions
        for i_session = 1:n_sessions
            if isfield(session_results(i_session), event_name) && ...
               isfield(session_results(i_session).(event_name), p_field)

                p_values = session_results(i_session).(event_name).(p_field);
                if isempty(p_values), continue; end

                if isempty(time_vector)
                    time_vector = session_results(i_session).(event_name).time_vector;
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
        y_max = max(y_max, max(get(h_axes(i_ax), 'YLim')));
    end
end
if y_max == 0, y_max = 1; end % Default if no data

linkaxes(h_axes(~strcmp(effect_terms(:,1),'---'), :), 'x');
set(h_axes(isgraphics(h_axes)), 'YLim', [0, y_max]);

for i_row = 1:n_rows
    if strcmp(effect_terms{i_row, 1}, '---'), continue; end
    for i_col = 1:n_cols
        ax = h_axes(i_row, i_col);
        if ~isgraphics(ax), continue; end
        if i_row == 1
            title(ax, strrep(analysis_plan.alignment_events{i_col}, '_', ' '));
        end
        if i_col == 1
            ylabel(ax, effect_terms{i_row, 1});
        else
            set(ax, 'YTickLabel', []);
        end
        if i_row < n_rows && ~strcmp(effect_terms{i_row+1,1}, '---')
            set(ax, 'XTickLabel', []);
        else
            xlabel(ax, 'Time (s)');
        end
    end
end

sgtitle(sprintf('Proportion of Significant Neurons (p < 0.05) for %s', brain_area_name), 'FontWeight', 'bold');

%% Save Figure
figures_dir = fullfile(project_root, 'figures', 'anova');
if ~exist(figures_dir, 'dir'), mkdir(figures_dir); end
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_anova_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end