%% plot_aggregated_decoding.m
%
% Generates a comprehensive summary figure visualizing both cross-factor and
% cross-time population decoding results from aggregated data.
%
% The figure is a multi-row grid of scatter plots where each point
% represents a single session. Different markers are used for each monkey,
% and error bars show the 95% confidence intervals for decoding accuracy.
%
% INPUTS:
%   aggregated_data - A struct array where each element corresponds to a
%                     session's aggregated analysis results. The structure
%                     is expected to contain decoding accuracy data.
%   brain_area_name - A string with the name of the brain area (e.g., 'SC')
%                     used for figure titles and filenames.
%   analysis_plan   - (Currently unused) A struct defining the analysis
%                     parameters. Included for future compatibility.
%
% ASSUMED DATA STRUCTURE:
% The function assumes that 'aggregated_data' has the following structure:
% aggregated_data(i_session).decoding_results
%   .Visual
%       .Reward
%           .accuracy       (scalar)
%           .accuracy_ci    ([lower, upper])
%           .generalization
%               .Salience
%                   .accuracy
%                   .accuracy_ci
%           .generalization_time
%               .Saccade
%                   .accuracy
%                   .accuracy_ci
%   .Saccade
%       ... (similar structure for Saccade epoch)
%
% Author: Jules
% Date: 2025-09-19
%

function plot_aggregated_decoding(aggregated_data, brain_area_name, analysis_plan)
%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Setup Figure
% Create a new figure window for the plots, making it large to accommodate
% the entire grid of subplots.
figure('Position', [100, 100, 1600, 1000], 'Name', ...
    sprintf('Aggregated Decoding Generalization: %s', brain_area_name));
hold on;

% Dynamically identify monkeys from session_ids. The monkey initial is
% assumed to be the first character of the session_id string.
session_ids = {aggregated_data.session_id};
monkey_ids = cellfun(@(x) x(1), session_ids, 'UniformOutput', false);
[unique_monkeys, ~, monkey_idx] = unique(monkey_ids);

% Define markers for each monkey. This can be extended if needed.
monkey_markers = {'o', 's', '^', 'd'};
if length(unique_monkeys) > length(monkey_markers)
    warning('More monkeys than defined markers; some will share markers.');
    marker_indices = mod(0:length(unique_monkeys)-1, length(monkey_markers)) + 1;
else
    marker_indices = 1:length(unique_monkeys);
end

%% Cross-Factor Generalization Plots
% This section creates a 2x6 grid of scatter plots (Rows 1-2 of the figure)
% to show generalization of decoders across different factors (e.g.,
% Reward, Salience) within the same epoch.

epochs = {'Visual', 'Saccade'};
factors = {'Reward', 'Salience', 'Identity', 'Probability'};
factor_pairs = nchoosek(factors, 2);
n_factor_pairs = size(factor_pairs, 1); % Should be 6

% Store handles for later formatting
h_axes_cross_factor = gobjects(length(epochs), n_factor_pairs);

for i_epoch = 1:length(epochs)
    epoch_name = epochs{i_epoch};

    for i_pair = 1:n_factor_pairs
        factor1 = factor_pairs{i_pair, 1};
        factor2 = factor_pairs{i_pair, 2};

        % Subplot position in a 4-row grid.
        plot_idx = (i_epoch - 1) * n_factor_pairs + i_pair;
        h_ax = mySubPlot([4, n_factor_pairs, plot_idx]);
        h_axes_cross_factor(i_epoch, i_pair) = h_ax;
        hold on;

        for i_session = 1:length(aggregated_data)
            % Check for existence of data before plotting
            if isfield(aggregated_data(i_session), 'decoding_results') && ...
               isfield(aggregated_data(i_session).decoding_results, epoch_name) && ...
               isfield(aggregated_data(i_session).decoding_results.(epoch_name), factor1) && ...
               isfield(aggregated_data(i_session).decoding_results.(epoch_name).(factor1), 'generalization') && ...
               isfield(aggregated_data(i_session).decoding_results.(epoch_name).(factor1).generalization, factor2)

                data = aggregated_data(i_session).decoding_results;

                % X-axis: Standard accuracy; Y-axis: Generalization accuracy
                x_val = data.(epoch_name).(factor1).accuracy;
                x_ci  = data.(epoch_name).(factor1).accuracy_ci;
                y_val = data.(epoch_name).(factor1).generalization.(factor2).accuracy;
                y_ci  = data.(epoch_name).(factor1).generalization.(factor2).accuracy_ci;

                marker_style = monkey_markers{marker_indices(monkey_idx(i_session))};

                % Plot point with X and Y error bars
                errorbar(x_val, y_val, ...
                    y_val - y_ci(1), y_ci(2) - y_val, ... % y error
                    x_val - x_ci(1), x_ci(2) - x_val, ... % x error
                    marker_style, 'MarkerSize', 6, 'LineWidth', 1, 'CapSize', 0);
            end
        end

        % Formatting for each subplot
        axis([0 1 0 1]);
        plot([0 1], [0 1], 'k:', 'LineWidth', 0.5); % Diagonal line
        axis square;
        title(sprintf('Train %s, Test %s', factor1, factor2));
    end
end

%% Cross-Time Generalization Plots
% This section creates a 2x4 grid of scatter plots (Rows 3-4 of the figure)
% to show generalization of decoders across time (from the visual to
% saccade epoch and vice versa).

% Row 3: Train Visual, Test Saccade
% Row 4: Train Saccade, Test Visual
train_epochs = {'Visual', 'Saccade'};
test_epochs = {'Saccade', 'Visual'};
n_factors = length(factors); % 'factors' is defined in the previous section

% Store handles for later formatting
h_axes_cross_time = gobjects(length(train_epochs), n_factors);

for i_dir = 1:length(train_epochs)
    train_epoch = train_epochs{i_dir};
    test_epoch = test_epochs{i_dir};

    for i_factor = 1:n_factors
        factor_name = factors{i_factor};

        % Position subplot in the 4x6 grid (rows 3 and 4, cols 1-4)
        plot_idx = (i_dir + 1) * n_factor_pairs + i_factor;
        h_ax = mySubPlot([4, n_factor_pairs, plot_idx]);
        h_axes_cross_time(i_dir, i_factor) = h_ax;
        hold on;

        for i_session = 1:length(aggregated_data)
            % Check for existence of data before plotting
            if isfield(aggregated_data(i_session), 'decoding_results') && ...
               isfield(aggregated_data(i_session).decoding_results, train_epoch) && ...
               isfield(aggregated_data(i_session).decoding_results.(train_epoch), factor_name) && ...
               isfield(aggregated_data(i_session).decoding_results.(train_epoch).(factor_name), 'generalization_time') && ...
               isfield(aggregated_data(i_session).decoding_results.(train_epoch).(factor_name).generalization_time, test_epoch)

                data = aggregated_data(i_session).decoding_results;

                % X-axis: Standard acc; Y-axis: Cross-time generalization acc
                x_val = data.(train_epoch).(factor_name).accuracy;
                x_ci  = data.(train_epoch).(factor_name).accuracy_ci;
                y_val = data.(train_epoch).(factor_name).generalization_time.(test_epoch).accuracy;
                y_ci  = data.(train_epoch).(factor_name).generalization_time.(test_epoch).accuracy_ci;

                marker_style = monkey_markers{marker_indices(monkey_idx(i_session))};

                errorbar(x_val, y_val, ...
                    y_val - y_ci(1), y_ci(2) - y_val, ...
                    x_val - x_ci(1), x_ci(2) - x_val, ...
                    marker_style, 'MarkerSize', 6, 'LineWidth', 1, 'CapSize', 0);
            end
        end

        % Formatting for each subplot
        axis([0 1 0 1]);
        plot([0 1], [0 1], 'k:', 'LineWidth', 0.5);
        axis square;
        title(sprintf('%s Decoder', factor_name));
    end
end

%% Finalize and Save Figure
% This section cleans up the figure by de-cluttering axes, adding shared
% labels, a legend, a main title, and then saves the final output as a PDF.

% --- De-clutter and Organize Axes ---
% Combine all axes handles into a single array for easier manipulation.
all_axes = gobjects(4, n_factor_pairs);
all_axes(1:2, :) = h_axes_cross_factor;
all_axes(3:4, 1:n_factors) = h_axes_cross_time;

% Hide unused subplot axes in the bottom two rows where there are no plots.
for r = 3:4
    for c = (n_factors + 1):n_factor_pairs
        if isgraphics(all_axes(r, c))
            set(all_axes(r, c), 'Visible', 'off');
        end
    end
end

valid_axes = all_axes(isgraphics(all_axes));
set(valid_axes, 'TickDir', 'out', 'Color', 'none', 'LineWidth', 1);

% Remove all tick labels initially, then restore them selectively to avoid
% clutter, following standard publication practices.
set(valid_axes, 'XTickLabel', [], 'YTickLabel', []);
set(all_axes(isgraphics(all_axes(:, 1))), 'YTickLabelMode', 'auto');

% Restore X-tick labels for the bottom-most plot in each column.
for c = 1:n_factor_pairs
    if c <= n_factors && isgraphics(all_axes(4, c))
        set(all_axes(4, c), 'XTickLabelMode', 'auto'); % Row 4 for cols 1-4
    elseif isgraphics(all_axes(2, c))
        set(all_axes(2, c), 'XTickLabelMode', 'auto'); % Row 2 for cols 5-6
    end
end

% --- Add Shared Labels, Row Labels, and Legend ---
% Create a "big" invisible axis to host shared labels.
big_ax = axes('Position', [0.05 0.05 0.9 0.9], 'Visible', 'off');
ylabel(big_ax, 'Generalization Accuracy', 'Visible', 'on', 'FontSize', 12);
xlabel(big_ax, 'Standard Accuracy', 'Visible', 'on', 'FontSize', 12);

% Add row labels using YLabels on the first plot of each row for clarity.
ylabel(h_axes_cross_factor(1,1), 'Visual Epoch');
ylabel(h_axes_cross_factor(2,1), 'Saccade Epoch');
ylabel(h_axes_cross_time(1,1), 'Train Vis -> Test Sac');
ylabel(h_axes_cross_time(2,1), 'Train Sac -> Test Vis');

% Create a legend for the monkey markers.
legend_handles = gobjects(1, length(unique_monkeys));
for i = 1:length(unique_monkeys)
    legend_handles(i) = plot(nan, nan, monkey_markers{marker_indices(i)}, ...
        'MarkerSize', 8, 'LineStyle', 'none', 'MarkerEdgeColor', 'k');
end
legend(legend_handles, unique_monkeys, 'Location', 'best', ...
    'Position', [0.8, 0.8, 0.1, 0.1], 'Box', 'off');

% --- Add Title and Save ---
sgtitle_str = sprintf('Cross-Factor and Cross-Time Decoding Generalization in %s', brain_area_name);
sgtitle(sgtitle_str, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');

% Define project root and figures directory using project conventions.
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
figures_dir = fullfile(project_root, 'figures');
if ~exist(figures_dir, 'dir'); mkdir(figures_dir); end

% Save the figure as a high-quality PDF.
fig = gcf;
fig_filename = fullfile(figures_dir, ...
    sprintf('aggregated_decoding_generalization_%s.pdf', brain_area_name));
pdfSave(fig_filename, fig.Position(3:4)/72, fig);

end
