%% plot_aggregated_baseline_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the baseline comparison analysis.
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-13
%

function plot_aggregated_baseline_comparison(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
% Define the project root, ensure the figures directory exists, and add the
% 'utils' directory to the path so that helper functions can be found.
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Figure and Plotting Setup
% Dynamically get the list of conditions from the data structure
if isfield(aggregated_sc_data, 'baseline_comparison')
    condition_names = fieldnames(aggregated_sc_data.baseline_comparison);
else
    condition_names = {};
end
n_conditions = length(condition_names);

if n_conditions == 0
    disp('No baseline_comparison data found to plot.');
    return;
end

fig = figure('Position', [100, 100, 400 * n_conditions, 600], ...
    'Color', 'w');
h_axes = gobjects(1, n_conditions);

% Define colors for SC and SNc
colors = richColors();
sc_color = colors(1,:);
snc_color = colors(2,:);
pos_alpha = 0.6; % Face alpha for positive modulation
neg_alpha = 0.3; % Face alpha for negative modulation

%% Main Plotting Loop
for i_cond = 1:n_conditions
    cond_name = condition_names{i_cond};

    % --- Data Validation ---
    if ~isfield(aggregated_sc_data, 'baseline_comparison') ...
            || ~isfield(aggregated_sc_data.baseline_comparison, cond_name)
        warning('plot_aggregated_baseline_comparison:no_sc_data', ...
            'No data for %s in SC struct.', cond_name);
        continue;
    end
    if ~isfield(aggregated_snc_data, 'baseline_comparison') ...
            || ~isfield(aggregated_snc_data.baseline_comparison, cond_name)
        warning('plot_aggregated_baseline_comparison:no_snc_data', ...
            'No data for %s in SNc struct.', cond_name);
        continue;
    end

    % --- Time Vector Loading ---
    time_vector = aggregated_sc_data.baseline_comparison.(cond_name).time_vector;

    % --- Data Extraction and Proportion Calculation ---
    sig_sc = aggregated_sc_data.baseline_comparison.(cond_name).sig;
    n_total_sc = size(sig_sc, 1);
    prop_sc_increase = sum(sig_sc == 1, 1) / n_total_sc;
    prop_sc_decrease = -sum(sig_sc == -1, 1) / n_total_sc; % Negative for plotting

    sig_snc = aggregated_snc_data.baseline_comparison.(cond_name).sig;
    n_total_snc = size(sig_snc, 1);
    prop_snc_increase = sum(sig_snc == 1, 1) / n_total_snc;
    prop_snc_decrease = -sum(sig_snc == -1, 1) / n_total_snc; % Negative for plotting

    % --- Plotting (SC and SNc on same axes) ---
    h_axes(1, i_cond) = mySubPlot([1, n_conditions, i_cond]);
    hold on;

    % Plot SC proportions
    h_sc_inc = barStairsFill(time_vector, zeros(size(prop_sc_increase)), prop_sc_increase);
    set(h_sc_inc(1), 'FaceColor', sc_color, 'EdgeColor', 'none', 'FaceAlpha', pos_alpha);
    delete(h_sc_inc(2:3));

    h_sc_dec = barStairsFill(time_vector, zeros(size(prop_sc_decrease)), prop_sc_decrease);
    set(h_sc_dec(1), 'FaceColor', sc_color, 'EdgeColor', 'none', 'FaceAlpha', neg_alpha);
    delete(h_sc_dec(2:3));

    % Plot SNc proportions
    h_snc_inc = barStairsFill(time_vector, zeros(size(prop_snc_increase)), prop_snc_increase);
    set(h_snc_inc(1), 'FaceColor', snc_color, 'EdgeColor', 'none', 'FaceAlpha', pos_alpha);
    delete(h_snc_inc(2:3));

    h_snc_dec = barStairsFill(time_vector, zeros(size(prop_snc_decrease)), prop_snc_decrease);
    set(h_snc_dec(1), 'FaceColor', snc_color, 'EdgeColor', 'none', 'FaceAlpha', neg_alpha);
    delete(h_snc_dec(2:3));

    % --- Formatting ---
    title_str = get_plot_title(cond_name);
    title(h_axes(1, i_cond), title_str, 'Interpreter', 'none');
    xlim([time_vector(1), time_vector(end)]);
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
    xlabel('Time from Event Onset (s)');
end

%% Figure Cleanup and Final Touches
% De-clutter axes
if n_conditions > 1
    set(h_axes(1, 2:end), 'YTickLabel', []);
end

% Add axis labels to outer plots
ylabel(h_axes(1, 1), 'Proportion of Neurons');

% Create a legend
h_leg_sc = patch(NaN, NaN, sc_color, 'FaceAlpha', pos_alpha);
h_leg_snc = patch(NaN, NaN, snc_color, 'FaceAlpha', pos_alpha);
legend([h_leg_sc, h_leg_snc], {'SC', 'SNc'}, 'Location', 'northeast', 'Box', 'off');

% set axes properties
allAx = findall(fig, 'Type', 'Axes');
set(allAx, 'TickDir', 'Out', 'Color', 'none', ...
    'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

% Save figure:
fig_filename = fullfile(figures_dir, 'aggregated_baseline_comparison.pdf');
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end

%% Local Helper Function
function title_str = get_plot_title(cond_name)
    % GET_PLOT_TITLE - Generates a plot title for a given condition.
    title_str = sprintf('Modulation during: %s', strrep(cond_name, '_', ' '));
end
