%% plot_aggregated_neuron_metrics.m
%
% Loads aggregated_neuron_metrics.mat and generates two figures:
% 1. A grid of grand average PSTHs for SC and SNc neurons.
% 2. A scatter plot of aggregated neuron metrics.
%
% Author: Jules
% Date: 2025-10-01
%

function plot_aggregated_neuron_metrics()
    %% Setup Paths
    % Add the 'utils' directory to the path so that helper functions can be
    % found.
    [script_dir, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(script_dir, 'utils'));

    %% Load Data
    % Define the path to the aggregated data file
    project_root = fileparts(script_dir);
    data_path = fullfile(project_root, 'data', 'processed', 'aggregated_neuron_metrics.mat');

    % Check if the file exists
    if ~exist(data_path, 'file')
        error('Aggregated data file not found: %s', data_path);
    end

    % Load the data
    fprintf('Loading aggregated data from %s...\n', data_path);
    load(data_path, 'psth_data', 'sc_metrics', 'snc_metrics', 'psth_time_vector');
    fprintf('Data loaded successfully.\n');

    %% Figure 1: Grand Average PSTHs
    fprintf('Generating Figure 1: Grand Average PSTHs...\n');

    % Define plot layout parameters
    brain_areas = {'sc', 'snc'};
    event_names = fieldnames(psth_data.sc);
    n_rows = numel(brain_areas);
    n_cols = numel(event_names);

    % Create figure and initialize axes array
    fig1 = figure('Position', [100, 100, 1500, 600]);
    ax1 = gobjects(n_rows, n_cols);

    % Loop through brain areas and events to create subplots
    for row = 1:n_rows
        area_name = brain_areas{row};

        for col = 1:n_cols
            event_name = event_names{col};

            % Create subplot using mySubPlot
            ax1(row, col) = mySubPlot([n_rows, n_cols, (row-1)*n_cols + col]);
            hold on;

            % Extract the corresponding PSTH data matrix
            % The data is stored in a struct with lowercase area names
            psth_matrix = psth_data.(area_name).(event_name);

            % Calculate the grand average PSTH (mean across neurons)
            mean_psth = mean(psth_matrix, 1);

            % Plot the grand average PSTH using the single-PSTH convention
            hBS = barStairsFill(psth_time_vector, zeros(size(mean_psth)), mean_psth);
            delete(hBS(2)); % Remove the baseline stair
            set(hBS(1), 'FaceColor', 'k', 'FaceAlpha', 0.8); % Set area fill
            set(hBS(3), 'Color', 'k'); % Set top line color

            title(event_name, 'Interpreter', 'none');
        end
    end

    % Formatting for Figure 1
    for row = 1:n_rows
        area_name_upper = upper(brain_areas{row});

        % Set common Y-axis limits for each row using outerLims
        outerLims(ax1(row, :), 'y');

        % Add Y-label only to the first plot in each row
        ylabel(ax1(row, 1), sprintf('%s Mean Firing Rate (Hz)', area_name_upper));

        % Remove Y-tick labels from all but the first column
        set(ax1(row, 2:end), 'YTickLabel', []);
    end

    for col = 1:n_cols
        event_name = event_names{col};

        % Add X-label only to the plots in the bottom row
        xlabel(ax1(n_rows, col), sprintf('Time from %s (s)', event_name));

        % Remove X-tick labels from all but the bottom row
        if n_rows > 1
            set(ax1(1:n_rows-1, col), 'XTickLabel', []);
        end
    end

    sgtitle('Figure 1: Grand Average PSTHs by Brain Area and Alignment Event', 'FontSize', 16);
    fprintf('Figure 1 generated.\n');

    %% Figure 2: Aggregated Neuron Metrics Scatter Plot
    fprintf('Generating Figure 2: Aggregated Neuron Metrics...\n');

    % Create figure
    fig2 = figure('Position', [100, 100, 1100, 500]);

    % Define plot layout parameters
    metrics_tables = {sc_metrics, snc_metrics};
    area_titles = {'Superior Colliculus', 'Substantia Nigra'};
    selected_colors = {'r', 'b'};

    % Store handles for legend
    legend_handles = gobjects(3, 1);
    legend_labels = {'Unclassified', 'Task-Modulated SC', 'Putative DA'};
    ax2 = gobjects(1, 2);

    for i = 1:2
        ax2(i) = mySubPlot([1, 2, i]);
        hold on;

        metrics_table = metrics_tables{i};

        % Plot Unselected Neurons
        unselected_idx = ~metrics_table.IsSelected;
        h_unselected = plot(metrics_table.WaveformDuration(unselected_idx), ...
             metrics_table.BaselineFR(unselected_idx), ...
             'o', 'MarkerSize', 4, 'MarkerEdgeColor', 'none', ...
             'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerFaceAlpha', 0.5);

        % Plot Selected Neurons
        selected_idx = metrics_table.IsSelected;
        h_selected = plot(metrics_table.WaveformDuration(selected_idx), ...
               metrics_table.BaselineFR(selected_idx), ...
               'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
               'MarkerFaceColor', selected_colors{i});

        % Store handles for legend
        if i == 1 % SC
            legend_handles(2) = h_selected;
        else % SNc
            legend_handles(3) = h_selected;
        end
    end
    legend_handles(1) = h_unselected; % From the last iteration, style is same

    % Formatting for Figure 2
    for i = 1:2
        set(ax2(i), 'YScale', 'log');
        title(ax2(i), area_titles{i});
        xlabel(ax2(i), 'Waveform Duration (ms)');
    end

    % Y-label only on the left plot
    ylabel(ax2(1), 'Baseline Firing Rate (Hz)');
    set(ax2(2), 'YTickLabel', []);

    % Add vertical reference line to SNc plot
    line(ax2(2), [0.6 0.6], get(ax2(2), 'YLim'), 'Color', 'k', 'LineStyle', '--');

    % Create a single legend for the figure
    legend(ax2(2), legend_handles, legend_labels, 'Location', 'northeastoutside');

    sgtitle('Figure 2: Waveform Duration vs. Baseline Firing Rate', 'FontSize', 16);
    fprintf('Figure 2 generated.\n');

    %% Save Figures
    fprintf('Saving figures...\n');
    figures_dir = fullfile(project_root, 'figures');
    if ~exist(figures_dir, 'dir')
        mkdir(figures_dir);
    end

    % Save Figure 1
    fig1_filename = fullfile(figures_dir, 'grand_average_psths.pdf');
    pdfSave(fig1, fig1_filename);
    fprintf('Saved Figure 1 to %s\n', fig1_filename);

    % Save Figure 2
    fig2_filename = fullfile(figures_dir, 'neuron_metrics_scatter.pdf');
    pdfSave(fig2, fig2_filename);
    fprintf('Saved Figure 2 to %s\n', fig2_filename);

    fprintf('All figures saved successfully.\n');

end