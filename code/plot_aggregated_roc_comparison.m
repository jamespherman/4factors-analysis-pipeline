%% plot_aggregated_roc_comparison.m
%
%   Generates a publication-quality summary figure that directly compares
%   SC and SNc population results from the bin-by-bin ROC comparison analysis.
%
% INPUTS:
%   aggregated_sc_data  - A struct containing aggregated data for SC.
%   aggregated_snc_data - A struct containing aggregated data for SNc.
%
% Author: Jules
% Date: 2025-09-12
%

function plot_aggregated_roc_comparison(aggregated_sc_data, aggregated_snc_data)

%% Setup Paths
% Define the project root, ensure the figures directory exists, and add the
% 'utils' directory to the path so that helper functions can be found.
project_root = fullfile(findOneDrive, 'Code', 'tokens-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
addpath(fullfile(project_root, 'code', 'utils'));

%% Figure and Plotting Setup
% Dynamically get the list of comparisons from the data structure
comparison_names = fieldnames(aggregated_sc_data.roc_comparison);
n_comparisons = length(comparison_names);

fig = figure('Position', [100, 100, 400 * n_comparisons, 800], ...
    'Color', 'w');
h_axes = gobjects(2, n_comparisons);

% Define colors for the two preferences (a > b and a < b)
colors = richColors;
posColor = colors(7,:);
negColor = colors(10,:);

%% Main Plotting Loop
for i_comp = 1:n_comparisons
    comp_name = comparison_names{i_comp};

    % --- Data Validation ---
    if ~isfield(aggregated_sc_data, 'roc_comparison') ...
            || ~isfield(aggregated_sc_data.roc_comparison, comp_name)
        warning('plot_aggregated_roc_comparison:no_sc_data', ...
            'No data for %s in SC struct.', comp_name);
        continue;
    end
    if ~isfield(aggregated_snc_data, 'roc_comparison') ...
            || ~isfield(aggregated_snc_data.roc_comparison, comp_name)
        warning('plot_aggregated_roc_comparison:no_snc_data', ...
            'No data for %s in SNc struct.', comp_name);
        continue;
    end

    % --- Time Vector Loading ---
    % The time vector is now self-contained in the aggregated data structure.
    time_vector = aggregated_sc_data.roc_comparison.(comp_name).time_vector;

    % --- Data Extraction and Count Calculation ---
    sig_sc = aggregated_sc_data.roc_comparison.(comp_name).sig;
    n_total_sc = size(sig_sc, 1);
    prop_sc_cond2 = sum(sig_sc == 1, 1);
    prop_sc_cond1 = -sum(sig_sc == -1, 1);

    sig_snc = aggregated_snc_data.roc_comparison.(comp_name).sig;
    n_total_snc = size(sig_snc, 1);
    prop_snc_cond2 = sum(sig_snc == 1, 1);
    prop_snc_cond1 = -sum(sig_snc == -1, 1);

    % --- Plotting SC Data (Top Row) ---
    h_axes(1, i_comp) = mySubPlot([2, n_comparisons, i_comp]);
    hold on;

    % Plot SC proportions: positive-going histogram first
    h_sc1 = barStairsFill(time_vector, zeros(size(prop_sc_cond2)), prop_sc_cond2);
    delete(h_sc1(2)); % Delete baseline
    set(h_sc1(1), 'FaceColor', posColor, 'EdgeColor', 'none');
    set(h_sc1(3), 'Color', posColor);

    % negative going next
    h_sc2 = barStairsFill(time_vector, zeros(size(prop_sc_cond1)), prop_sc_cond1);
    delete(h_sc2(2)); % Delete baseline
    set(h_sc2(1), 'FaceColor', negColor, 'EdgeColor', 'none');
    set(h_sc2(3), 'Color', negColor);

    % --- Get Plot Labels ---
    % We can use either aggregated_sc_data or aggregated_snc_data, as the
    % alignment event should be the same for a given comparison.
    [title_str, xlabel_str] = ...
        get_plot_labels(comp_name, aggregated_sc_data);

    % Formatting for SC plots
    title(h_axes(1, i_comp), title_str, 'Interpreter', 'none');
    xlim([time_vector(1), time_vector(end)]);
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');

    % --- Plotting SNc Data (Bottom Row) ---
    h_axes(2, i_comp) = mySubPlot([2, n_comparisons, i_comp + ...
        n_comparisons]);
    hold on;

    % Plot SNc proportions, positive first
    h_snc1 = barStairsFill(time_vector, zeros(size(prop_snc_cond2)), ...
        prop_snc_cond2);
    delete(h_snc1(2));
    set(h_snc1(1), 'FaceColor', posColor, 'EdgeColor', 'none');
    set(h_snc1(3), 'Color', posColor);

    % negaitve next:
    h_snc2 = barStairsFill(time_vector, zeros(size(prop_snc_cond1)), ...
        prop_snc_cond1);
    delete(h_snc2(2));
    set(h_snc2(1), 'FaceColor', negColor, 'EdgeColor', 'none');
    set(h_snc2(3), 'Color', negColor);

    % Formatting for SNc plots
    xlabel(h_axes(2, i_comp), xlabel_str);
    xlim([time_vector(1), time_vector(end)]);
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', '--');
    line([0, 0], ylim, 'Color', 'k', 'LineStyle', '--');
end

%% Figure Cleanup and Final Touches
% De-clutter axes per AGENTS.md instructions
 % Remove x-labels from top row
set(h_axes(1, :), 'XTickLabel', []);
if n_comparisons > 1
     % Remove y-labels from all but the first column
    set(h_axes(:, 2:end), 'YTickLabel', []);
end

% Add axis labels to outer plots
ylabel(h_axes(1, 1), 'Count of Neurons (SC)');
ylabel(h_axes(2, 1), 'Count of Neurons (SNc)');

% set axes properties
allAx = findall(fig, 'Type', 'Axes');

% set common y-limits:
[~, yLims] = outerLims(allAx);
set(allAx, 'YLim', yLims, 'TickDir', 'Out', 'Color', 'none', ...
    'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

% Save figure:
fig_filename = fullfile(figures_dir, 'aggregated_roc_comparison.pdf');
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end

%% Local Helper Functions
function [title_str, xlabel_str] = get_plot_labels(comp_name, aggregated_data)
    % GET_PLOT_LABELS - Generates plot titles and labels for a given comparison.
    %
    % INPUTS:
    %   comp_name       - The name of the comparison (e.g., 'RPE_at_Outcome').
    %   aggregated_data - The aggregated data struct for one population.
    %
    % OUTPUTS:
    %   title_str       - A formatted string for the plot title.
    %   xlabel_str      - A formatted string for the x-axis label.

    % Generate title string based on comparison name
    switch comp_name
        case 'Dist_at_Cue'
            title_str = 'Preference: Normal vs. Uniform at Cue Onset';
        case 'RPE_at_Outcome'
            title_str = 'Preference: Common vs. Rare High Reward at Outcome';
        case 'RPE_at_Reward'
            title_str = 'Preference: Common vs. Rare High Reward at Reward';
        case 'SPE_at_Outcome'
            title_str = 'Preference: Surprise vs. No Surprise at Outcome';
        otherwise
            % Default title format if no specific case matches
            title_str = strrep(comp_name, '_', ' ');
    end

    % Define the regular expression pattern
    pattern = '^[^_]*_[^_]*_(.*)$';
    matches = regexp(comp_name, pattern, 'tokens');


    % Generate xlabel string using the alignment event
    align_event = matches{1}{1};
    xlabel_str = sprintf('Time from %s Onset (s)', strrep(align_event, '_', ' '));
end
