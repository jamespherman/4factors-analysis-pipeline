function generate_neuron_summary_pdf(session_data, ...
    selected_neurons, unique_id, output_dir)
% GENERATE_NEURON_SUMMARY_PDF Generates a multi-page PDF with
% diagnostic plots for each neuron.
%
%   GENERATE_NEURON_SUMMARY_PDF(session_data, selected_neurons, ...
%   unique_id)
%   generates a multi-page PDF file where each page contains a summary
%   of a single neuron, conforming to the standardized `session_data`
%   structure.
%
%   Inputs:
%       session_data     - A struct for a single session, containing
%                          fields like 'spikes' and 'cluster_info'.
%       selected_neurons - A logical vector (nClusters x 1) indicating
%                          which neurons were selected by a screening
%                          process.
%       unique_id        - A string used to name the output PDF file.
%       output_dir       - A string used to indicate where the output PDF
%                          file will be stored.
%
%   Output:
%       The function saves a PDF file named
%       '[unique_id]_neuron_diagnostics.pdf' in the 'figures/output_dir/'
%       directory.

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can
% be found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Analysis Plan
% Retrieve the analysis plan and session-specific condition masks.
analysis_plan = define_task_conditions();
condition_masks = define_task_conditions(session_data);

%% Constants
% Define constants for plot layout and analysis parameters.
N_ROWS = 12;
N_COLS = 6;
PSTH_WINDOW = [-0.5, 1.0]; % 1.5s window for PSTHs
SLIDING_BIN_WIDTH = 0.2;   % 200ms sliding bin
SLIDING_STEP_SIZE = 0.1;   % 100ms step size

%% Layout Configuration
% Define column counts for different sections of the plot
N_COLS1 = 6; % For general diagnostics
N_COLS2 = 4; % For PSTH/raster pairs

% Pre-calculate subplot indices for a complex layout
% This allows mixing plots from different conceptual grids on the same figure
grid1_indices = reshape(1:(N_COLS1 * N_ROWS), N_COLS1, N_ROWS)';
grid2_indices = reshape(1:(N_COLS2 * N_ROWS), N_COLS2, N_ROWS)';

% Indices for the left-side plots (Waveform, ISI) using the 6-column grid
left_plot_indices = grid1_indices(:, 1:2);

% Indices for the right-side plots (PSTHs) using the 4-column grid
right_plot_indices = grid2_indices(:, 3:4);

%% Configuration
% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

%% --- Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info.cluster_id);
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

% This space is intentionally left blank.

% Loop through each cluster to create a page of plots
for i_cluster = 1:nClusters

    % Create a new figure for the PDF
    fig = figure('Color', 'w', 'MenuBar', 'None', 'ToolBar', 'None');

    % --- Layout Configuration ---
    top_panel_grid = [N_ROWS, N_COLS];
    psth_grid      = [N_ROWS, 4];


    % define CLUSTER ID:
    cluster_id = cluster_ids(i_cluster);
    
    % define output file name:
    output_filename = fullfile(output_dir, ...
        [unique_id '_' sprintf('%03d',i_cluster) ...
        '_neuron_diagnostics.pdf']);

    % Add a title for the page
    sgtitle(sprintf('Neuron Diagnostic Summary: Cluster %d', cluster_id));

    spike_times = all_spike_times(all_spike_clusters == cluster_id);

    %% --- Top Panel Diagnostics ---
    % This section is now outside the conditional `do_spatial_tuning_plot`
    % block to avoid code duplication. It uses the pre-calculated indices
    % to place plots in a complex, multi-grid layout.

    % Waveform Plot (Position 1)
    ax_wf = mySubPlot([N_ROWS, N_COLS1, left_plot_indices(1,1)]);
    set(ax_wf, 'Tag', 'Waveform_Axis');
    if isfield(session_data.spikes, 'wfMeans') && ...
            numel(session_data.spikes.wfMeans) >= i_cluster
        plot(ax_wf, session_data.spikes.wfMeans{i_cluster}');
        title(ax_wf, 'Mean Waveform');
        xlabel(ax_wf, 'Samples');
        ylabel(ax_wf, 'Amplitude (uV)');
        axis(ax_wf, 'tight');
    else
        text(0.5, 0.5, 'No waveform', 'Parent', ax_wf, ...
            'HorizontalAlignment', 'center');
    end

    % ISI Histogram (Position 2)
    ax_isi = mySubPlot([N_ROWS, N_COLS1, left_plot_indices(1,2)]);
    set(ax_isi, 'Tag', 'ISI_Axis');
    if numel(spike_times) > 1
        isi = diff(spike_times) * 1000; % ms
        histogram(ax_isi, isi, 'EdgeColor', 'k', ...
            'FaceColor', [0.5 0.5 0.5]);
        set(ax_isi, 'XScale', 'log');
        title(ax_isi, 'ISI Histogram');
        xlabel(ax_isi, 'ISI (ms, log)');
        ylabel(ax_isi, 'Count');
    else
        text(0.5, 0.5, 'No ISI', 'Parent', ax_isi, ...
            'HorizontalAlignment', 'center');
    end

    %% --- Dynamic PSTH and Raster Plots ---
    n_diag_plots = numel(analysis_plan.diagnostic_plots);

    for i_plot = 1:n_diag_plots
        plot_def = analysis_plan.diagnostic_plots(i_plot);
        event_name = plot_def.event;
        plot_title = plot_def.title;
        conditions_to_compare = plot_def.conditions_to_compare;

        % Create axes for the raster and PSTH plots
        ax_raster = mySubPlot([N_ROWS, N_COLS2, right_plot_indices(1, i_plot)]);
        set(ax_raster, 'Tag', 'Raster_Axis');
        ax_psth = mySubPlot([N_ROWS, N_COLS2, right_plot_indices(1, i_plot) + N_COLS2]);
        set(ax_psth, 'Tag', 'PSTH_Axis');

        hold(ax_psth, 'on');

        colors = richColors(numel(conditions_to_compare));
        all_event_times_for_raster = [];

        for i_cond = 1:numel(conditions_to_compare)
            cond_name = conditions_to_compare{i_cond};

            % Get trial indices for the current condition
            trial_mask = condition_masks.(cond_name);

            % Get valid event times for these trials
            event_times = session_data.eventTimes.(event_name)(trial_mask);
            event_times = event_times(~isnan(event_times) & event_times > 0);

            if isempty(event_times)
                continue;
            end

            all_event_times_for_raster = [all_event_times_for_raster; event_times];

            % Calculate PSTH for the current condition
            [psth, bin_centers, ~] = alignAndBinSpikes(spike_times, ...
                event_times, PSTH_WINDOW(1), PSTH_WINDOW(2), ...
                SLIDING_BIN_WIDTH, SLIDING_STEP_SIZE);

            n_trials = size(psth, 1);
            if n_trials == 0
                continue;
            end

            mean_rate = sum(psth, 1, 'omitnan') / (n_trials * SLIDING_BIN_WIDTH);

            % Plot the mean PSTH for the current condition
            hBS = barStairsFill(bin_centers, zeros(size(mean_rate)), mean_rate);
            delete(hBS(1:2)); % Keep only the line for overlay
            set(hBS(3), 'Color', colors(i_cond,:));
        end

        % Create a single raster plot for all conditions
        if ~isempty(all_event_times_for_raster)
            [psth_raster, bin_centers_raster, ~] = alignAndBinSpikes(spike_times, ...
                all_event_times_for_raster, PSTH_WINDOW(1), PSTH_WINDOW(2), ...
                SLIDING_BIN_WIDTH, SLIDING_STEP_SIZE);

            imagesc(ax_raster, bin_centers_raster, 1:size(psth_raster, 1), psth_raster);
            colormap(ax_raster, flipud(bone(64)));
            set(ax_raster, 'YDir', 'normal');
            title(ax_raster, plot_title);
            set(ax_raster, 'xticklabel', {[]});
        else
            text(0.5, 0.5, 'No data for raster', 'Parent', ax_raster, 'HorizontalAlignment', 'center');
        end

        % Finalize PSTH plot appearance
        xlabel(ax_psth, sprintf('Time from %s (s)', plot_title));
        if i_plot == 1
            ylabel(ax_psth, 'Firing Rate (Hz)');
        else
            set(ax_psth, 'yticklabel', {[]});
        end
        xline(ax_psth, 0, 'r--');

        % Link axes and set limits
        linkaxes([ax_raster, ax_psth], 'x');
        xlim(ax_raster, PSTH_WINDOW);
    end

    %% --- Summary Information ---
    ax_summary = mySubPlot([top_panel_grid, N_COLS * (N_ROWS - 1) + 1]);
    set(ax_summary, 'Tag', 'Summary_Axis');
    axis(ax_summary, 'off')

    screening_status = 'Not Selected';
    if numel(selected_neurons) >= i_cluster && selected_neurons(i_cluster)
        screening_status = 'Selected';
    end

    phy_quality = 'N/A';
    info_row = cluster_info.cluster_id == cluster_id;
    if any(info_row) && isfield(cluster_info, 'group')
        phy_quality = cluster_info.group(info_row);
    end

    baseline_fr = session_data.metrics.baseline_frs(i_cluster);
    wf_metrics = session_data.metrics.wf_metrics(i_cluster);
    wf_duration = sprintf('%.2f ms', wf_metrics.peak_trough_ms);

    summary_text = {
        sprintf('Cluster ID: %d', cluster_id), ...
        sprintf('Phy Quality: %s', phy_quality), ...
        sprintf('Screening Status: %s', screening_status), ...
        sprintf('Baseline FR: %.2f Hz', baseline_fr), ...
        sprintf('Waveform Duration: %s', wf_duration)
    };

    text(0.1, 0.5, summary_text, 'Parent', ax_summary, ...
        'VerticalAlignment', 'middle', 'FontSize', 10);
    titleObj = title(ax_summary, 'Summary Information');
    set(titleObj, 'Position', [0.5 1.5 0.5])

    %% --- Final Formatting ---
    % Find all PSTH axes by tag
    psth_axes = findobj(fig, 'Tag', 'PSTH_Axis');

    % Get 'outer limits' for all PSTH axes
    [~, yLims] = outerLims(psth_axes);

    % General appearance modifications
    set(findall(fig, 'Type', 'Axes'), 'XColor', 'k', 'YColor', 'k', ...
        'TickDir', 'Out', 'LineWidth', 1, 'Color', 'w', 'Box', 'Off');
    set(findall(fig, 'Type', 'Text'), 'Color', 0.2*[1 1 1], 'FontSize', 8);
    set(fig, 'Position', [150 50 750 1000]);

    % Synchronize Y-limits for all PSTH axes
    set(psth_axes, 'YLim', yLims);

    drawnow;

    % Appending only works with postscript files in MATLAB so we will
    % create one PDF per unit:
    pdfSave(output_filename, fig.Position(3:4)/72, fig);

    % give feedback:
    giveFeed(['Wrote PDF for unit ' sprintf('%03d',i_cluster)])

    % Close the figure
    close(fig);
end

fprintf('Diagnostic PDF saved to %s\n', output_filename);

end