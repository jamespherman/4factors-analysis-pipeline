function generate_neuron_summary_pdf(session_data, ...
    selected_neurons, unique_id, output_dir, varargin)
% GENERATE_NEURON_SUMMARY_PDF Generates single-page PDF summaries.
%
%   GENERATE_NEURON_SUMMARY_PDF(session_data, selected_neurons, unique_id,
%   output_dir) generates a single-page PDF file for each neuron.
%
%   Layout (3 independent column groups using viewport partitioning):
%       Group 1 (x: 0.02-0.08): Waveform + ISI + Text info (stacked)
%       Group 2 (x: 0.13-0.25): Spatial tuning (2 cols × 4 rows)
%       Group 3 (x: 0.30-0.98): Factor sensitivity (12 cols × 4 rows)
%
%   Each PSTH panel has a raster plot above it. Factor sensitivity
%   rasters/PSTHs use color-coded conditions (blue=High/Face,
%   orange=Low/Non-face).
%
%   Output:
%       PDF files named '{unique_id}_{cluster_id:03d}_summary.pdf'

%% Setup Paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Parse Optional Arguments
p = inputParser;
addParameter(p, 'ClusterIndex', [], @isnumeric);
parse(p, varargin{:});
cluster_index_to_process = p.Results.ClusterIndex;

%% Constants
N_ROWS = 4;
N_LOCATIONS = 4;

% PSTH parameters
PSTH_WINDOW = [-0.5, 1.0];
SLIDING_BIN_WIDTH = 0.025;
SLIDING_STEP_SIZE = 0.025;

% Colors from richColors palette
colors_palette = richColors('matrix');
COLOR_HIGH = colors_palette(6, :);   % Bright Blue [26, 133, 255]/255
COLOR_LOW = colors_palette(12, :);   % Dark Orange [153, 79, 0]/255
COLOR_BLACK = [0 0 0];

% Layout parameters for the 3 column groups (viewport partitioning)
% Group 1: Waveform/ISI/Text (single column)
G1_LEFT_MARGIN = 0.02;
G1_WIDTH = 0.06;

% Group 2: Spatial tuning (2 columns × 8 rows)
G2_LEFT_MARGIN = 0.13;
G2_WIDTH = 0.12;

% Group 3: Factor sensitivity (12 columns × 8 rows)
G3_LEFT_MARGIN = 0.30;
G3_WIDTH = 0.68;

% Common vertical parameters
GRID_HEIGHT = 0.88;
GRID_BOTTOM_MARGIN = 0.07;
HEIGHT_SPACING = 0.02;
G2_WIDTH_SPACING = 0.008;
G3_WIDTH_SPACING = 0.008;

% X-tick labels for waveform plots
xticklabels_wf = arrayfun(@num2str, round((0:20:80)./30, 2), ...
    'UniformOutput', false);

%% Setup and Data Extraction
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = numel(cluster_ids);
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

% Get condition masks for factor comparisons
codes = initCodes;
[conditions, ~] = define_task_conditions(session_data);

% Get 4factors trial indices (valid rewarded trials)
fourfactors_task_code = codes.uniqueTaskCode_gSac_4factors;
fourfactors_master_mask = session_data.trialInfo.taskCode == ...
    fourfactors_task_code & ~cellfun(@isempty, ...
    session_data.eventTimes.rewardCell);
fourfactors_trial_indices = find(fourfactors_master_mask);

%% Get target locations from 4factors trials
targetLocIdx_all = session_data.trialInfo.targetLocIdx;
targetTheta_all = session_data.trialInfo.targetTheta;

targetLocIdx_4f = targetLocIdx_all(fourfactors_trial_indices);
targetTheta_4f = targetTheta_all(fourfactors_trial_indices) / 10;

% Get unique locations sorted by targetLocIdx
[unique_locs, ~, loc_group_idx] = unique(targetLocIdx_4f);
n_locations = min(N_LOCATIONS, length(unique_locs));

% Get representative theta for each location
loc_thetas = zeros(n_locations, 1);
for i_loc = 1:n_locations
    loc_trials = loc_group_idx == i_loc;
    loc_thetas(i_loc) = median(targetTheta_4f(loc_trials));
end

%% Determine which clusters to process
if ~isempty(cluster_index_to_process)
    clusters_to_process = cluster_index_to_process;
else
    clusters_to_process = 1:nClusters;
end

%% Loop through each cluster
for i_cluster = clusters_to_process

    cluster_id = cluster_ids(i_cluster);
    spike_times = all_spike_times(all_spike_clusters == cluster_id);

    % Define output filename
    output_filename = fullfile(output_dir, ...
        sprintf('%s_%03d_summary.pdf', unique_id, cluster_id));

    % Create figure at 90% screen size
    screenSize = get(0, 'ScreenSize');
    figSize = [screenSize(3) * 0.9, screenSize(4) * 0.9];
    fig = figure('Color', 'w', 'MenuBar', 'None', 'ToolBar', 'None', ...
        'Position', [screenSize(3)*0.05, screenSize(4)*0.05, ...
        figSize(1), figSize(2)]);

    %% --- Header ---
    if numel(selected_neurons) >= i_cluster && selected_neurons(i_cluster)
        screening_status = 'SELECTED';
        status_color = [0.2 0.6 0.2];
    else
        screening_status = 'NOT SELECTED';
        status_color = [0.6 0.2 0.2];
    end

    % Status indicator in top-right
    annotation('textbox', [0.85 0.96 0.14 0.03], ...
        'String', screening_status, ...
        'FontSize', 10, 'FontWeight', 'bold', 'EdgeColor', status_color, ...
        'BackgroundColor', status_color, 'Color', 'w', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'LineWidth', 1.5);

    %% --- Group 1: Waveform + ISI + Text Info ---
    opts_g1 = {'LeftMargin', G1_LEFT_MARGIN, 'Width', G1_WIDTH, ...
        'Height', GRID_HEIGHT, 'BottomMargin', GRID_BOTTOM_MARGIN, ...
        'HeightSpacing', HEIGHT_SPACING};

    % Waveform: 3/4 of column height at top position
    ax_wf = mySubPlot([5/3, 1, 1], opts_g1{:});
    if isfield(session_data.spikes, 'wfMeans') && ...
            numel(session_data.spikes.wfMeans) >= i_cluster
        ow = offsetArrayWaveform(...
            session_data.spikes.wfMeans{i_cluster}, 0.025);
        or = range(ow(:));
        om = mean(ow(:));
        yLim = om + [-1 1] * (1.1 * or) / 2;
        % Plot all channels in black
        plot(ax_wf, ow', 'Color', COLOR_BLACK);
        title(ax_wf, 'Waveform', 'FontSize', 8);
        set(ax_wf, 'YLim', yLim, 'Box', 'Off', 'TickDir', 'Out', ...
            'YTickLabel', '', 'XLim', [0 83], 'XTick', 0:20:80, ...
            'XTickLabel', xticklabels_wf, 'FontSize', 7);
        xlabel(ax_wf, 'Time (ms)', 'FontSize', 7);
    else
        text(0.5, 0.5, 'No waveform', 'Parent', ax_wf, ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
        axis(ax_wf, 'off');
    end

    % ISI histogram: 1/8 of column height at row 7 of 8
    ax_isi = mySubPlot([8, 1, 7], opts_g1{:});
    if numel(spike_times) > 1
        isi = diff(spike_times) * 1000;
        histogram(ax_isi, isi, 30, 'EdgeColor', 'none', ...
            'FaceColor', [0.5 0.5 0.5]);
        set(ax_isi, 'XScale', 'log', 'FontSize', 6);
        title(ax_isi, 'ISI', 'FontSize', 7);
        refrac_violations = sum(isi < 1.5) / length(isi) * 100;
        text(0.95, 0.85, sprintf('%.1f%%<1.5ms', refrac_violations), ...
            'Parent', ax_isi, 'Units', 'normalized', ...
            'HorizontalAlignment', 'right', 'FontSize', 5);
    else
        text(0.5, 0.5, 'No ISI', 'Parent', ax_isi, ...
            'HorizontalAlignment', 'center', 'FontSize', 7);
        axis(ax_isi, 'off');
    end

    % Text info: 1/8 of column height at row 8 of 8 (bottom)
    ax_info = mySubPlot([8, 1, 8], opts_g1{:});
    axis(ax_info, 'off');

    % Get neuron metadata
    phy_quality = 'N/A';
    info_row = cluster_info.cluster_id == cluster_id;
    if any(info_row) && isfield(cluster_info, 'group')
        phy_quality = cluster_info.group(info_row);
        if iscell(phy_quality)
            phy_quality = phy_quality{1};
        end
    end

    baseline_fr = session_data.metrics.baseline_frs(i_cluster);
    wf_metrics = session_data.metrics.wf_metrics(i_cluster);
    wf_duration = wf_metrics.peak_trough_ms;
    brain_area = session_data.metadata.brain_area;

    info_str = sprintf(['%s\nCluster %d | %s\n' ...
        'FR:%.1f Hz | WF:%.2f ms\nPhy:%s | n4F:%d'], ...
        unique_id, cluster_id, brain_area, ...
        baseline_fr, wf_duration, phy_quality, ...
        length(fourfactors_trial_indices));

    text(0.5, 0.5, info_str, 'Parent', ax_info, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 6, 'Interpreter', 'none');

    %% --- Precompute all PSTH and Raster data ---
    % Event times extraction
    events = {'targetOn', 'saccadeOnset', 'reward'};
    event_times_all = cell(1, 3);
    for i_ev = 1:3
        if strcmp(events{i_ev}, 'reward')
            reward_cell = session_data.eventTimes.rewardCell;
            et = NaN(size(reward_cell));
            for i_t = 1:numel(reward_cell)
                if ~isempty(reward_cell{i_t})
                    et(i_t) = reward_cell{i_t}(1);
                end
            end
            event_times_all{i_ev} = et(fourfactors_trial_indices);
        else
            event_times_all{i_ev} = session_data.eventTimes.(events{i_ev})(...
                fourfactors_trial_indices);
        end
    end

    % Spatial tuning PSTHs (by location, all factors collapsed)
    spatial_psth_data = cell(n_locations, 2);
    for i_loc = 1:n_locations
        loc_mask = loc_group_idx == i_loc;
        for i_ev = 1:2
            event_times = event_times_all{i_ev};
            if sum(loc_mask) >= 3
                [psth, bin_centers, ~] = alignAndBinSpikes(spike_times, ...
                    event_times(loc_mask), PSTH_WINDOW(1), ...
                    PSTH_WINDOW(2), SLIDING_BIN_WIDTH, SLIDING_STEP_SIZE);
                rate = mean(psth, 1, 'omitnan') / SLIDING_BIN_WIDTH;
                spatial_psth_data{i_loc, i_ev} = struct(...
                    'rate', rate, 'bin_centers', bin_centers, ...
                    'n_trials', sum(loc_mask), 'loc_mask', loc_mask);
            end
        end
    end

    % Factor definitions
    factors = {
        'Reward', 'is_high_reward', 'is_low_reward', [];
        'Salience', 'is_high_salience', 'is_low_salience', ...
            'is_bullseye_target';
        'Probability', 'is_high_probability', 'is_low_probability', [];
        'Identity', 'is_face_target', 'is_nonface_target', 'is_image_target'
    };

    % Factor PSTHs and raster trial masks
    factor_psth_data = cell(4, n_locations, 3);
    factor_trial_masks = cell(4, n_locations);  % Store condition masks per factor/location

    for i_factor = 1:4
        cond1_mask = conditions.(factors{i_factor, 2});
        cond2_mask = conditions.(factors{i_factor, 3});
        if ~isempty(factors{i_factor, 4})
            filter_mask = conditions.(factors{i_factor, 4});
            cond1_mask = cond1_mask & filter_mask;
            cond2_mask = cond2_mask & filter_mask;
        end

        for i_loc = 1:n_locations
            loc_mask = loc_group_idx == i_loc;
            c1_loc = cond1_mask & loc_mask;
            c2_loc = cond2_mask & loc_mask;

            % Store trial masks for rasters
            factor_trial_masks{i_factor, i_loc} = struct(...
                'cond1', c1_loc, 'cond2', c2_loc);

            for i_ev = 1:3
                event_times = event_times_all{i_ev};
                n1 = sum(c1_loc);
                n2 = sum(c2_loc);

                if n1 >= 3 && n2 >= 3
                    [psth1, bin_centers, ~] = alignAndBinSpikes(...
                        spike_times, event_times(c1_loc), ...
                        PSTH_WINDOW(1), PSTH_WINDOW(2), ...
                        SLIDING_BIN_WIDTH, SLIDING_STEP_SIZE);
                    [psth2, ~, ~] = alignAndBinSpikes(...
                        spike_times, event_times(c2_loc), ...
                        PSTH_WINDOW(1), PSTH_WINDOW(2), ...
                        SLIDING_BIN_WIDTH, SLIDING_STEP_SIZE);

                    rate1 = mean(psth1, 1, 'omitnan') / SLIDING_BIN_WIDTH;
                    rate2 = mean(psth2, 1, 'omitnan') / SLIDING_BIN_WIDTH;

                    factor_psth_data{i_factor, i_loc, i_ev} = struct(...
                        'rate1', rate1, 'rate2', rate2, ...
                        'bin_centers', bin_centers, ...
                        'n1', n1, 'n2', n2, ...
                        'cond1_mask', c1_loc, 'cond2_mask', c2_loc);
                end
            end
        end
    end

    %% --- Compute shared Y-limits ---
    spatial_ymax = 0;
    for i_loc = 1:n_locations
        for i_ev = 1:2
            if ~isempty(spatial_psth_data{i_loc, i_ev})
                spatial_ymax = max(spatial_ymax, ...
                    max(spatial_psth_data{i_loc, i_ev}.rate));
            end
        end
    end
    spatial_ylim = [0, spatial_ymax * 1.1 + 1];

    factor_ymax = 0;
    for i_factor = 1:4
        for i_loc = 1:n_locations
            for i_ev = 1:3
                if ~isempty(factor_psth_data{i_factor, i_loc, i_ev})
                    d = factor_psth_data{i_factor, i_loc, i_ev};
                    factor_ymax = max([factor_ymax, max(d.rate1), ...
                        max(d.rate2)]);
                end
            end
        end
    end
    factor_ylim = [0, factor_ymax * 1.1 + 1];

    %% --- Group 2: Spatial Tuning (2 cols × 8 rows with rasters) ---
    % Using 8-row grid: odd rows for rasters, even rows for PSTHs
    % Each location spans 2 rows: row 2k-1 = raster, row 2k = PSTH

    event_labels = {'Target On', 'Saccade'};

    % Options for Group 2: 8 rows × 2 columns grid
    opts_g2 = {'LeftMargin', G2_LEFT_MARGIN, 'Width', G2_WIDTH, ...
        'Height', GRID_HEIGHT, 'BottomMargin', GRID_BOTTOM_MARGIN, ...
        'HeightSpacing', 0.008, 'WidthSpacing', G2_WIDTH_SPACING};

    for i_loc = 1:n_locations
        for i_ev = 1:2
            col = i_ev;

            % Calculate grid indices for 8-row × 2-column grid
            % Raster row: 1, 3, 5, 7 for locations 1, 2, 3, 4
            % PSTH row: 2, 4, 6, 8 for locations 1, 2, 3, 4
            raster_row = (i_loc - 1) * 2 + 1;
            psth_row = (i_loc - 1) * 2 + 2;

            % Calculate subplot indices (row-major order)
            raster_idx = (raster_row - 1) * 2 + col;
            psth_idx = (psth_row - 1) * 2 + col;

            % Create raster axes first (top of pair)
            ax_raster = mySubPlot([8, 2, raster_idx], opts_g2{:});

            % Create PSTH axes last so it's current for barStairsFill
            ax_psth = mySubPlot([8, 2, psth_idx], opts_g2{:});

            if ~isempty(spatial_psth_data{i_loc, i_ev})
                d = spatial_psth_data{i_loc, i_ev};

                % Plot raster (all trials in black)
                hold(ax_raster, 'on');
                loc_trials = find(d.loc_mask);
                for i_t = 1:length(loc_trials)
                    trial_idx = loc_trials(i_t);
                    aligned_spikes = spike_times - ...
                        event_times_all{i_ev}(trial_idx);
                    spikes_in_window = aligned_spikes(...
                        aligned_spikes >= PSTH_WINDOW(1) & ...
                        aligned_spikes <= PSTH_WINDOW(2));
                    if ~isempty(spikes_in_window)
                        plot(ax_raster, spikes_in_window, ...
                            i_t * ones(size(spikes_in_window)), ...
                            '.', 'Color', COLOR_BLACK, 'MarkerSize', 1);
                    end
                end
                xline(ax_raster, 0, 'k--', 'LineWidth', 0.5);
                set(ax_raster, 'XLim', PSTH_WINDOW, ...
                    'YLim', [0 length(loc_trials)+1], ...
                    'XTickLabel', [], 'YTick', [], ...
                    'TickDir', 'Out', 'Box', 'Off');

                % Plot PSTH
                hold(ax_psth, 'on');
                hBS = barStairsFill(d.bin_centers, ...
                    zeros(size(d.rate)), d.rate);
                delete(hBS(1:2));
                set(hBS(3), 'Color', COLOR_BLACK, 'LineWidth', 1);
                xline(ax_psth, 0, 'k--', 'LineWidth', 0.5);
                set(ax_psth, 'XLim', PSTH_WINDOW, 'YLim', spatial_ylim, ...
                    'TickDir', 'Out', 'Box', 'Off', 'FontSize', 6);

                % Title (top row only)
                if i_loc == 1
                    title(ax_raster, event_labels{i_ev}, 'FontSize', 7);
                end
            else
                axis(ax_raster, 'off');
                text(0.5, 0.5, 'No data', 'Parent', ax_psth, ...
                    'HorizontalAlignment', 'center', 'FontSize', 6);
                set(ax_psth, 'XLim', PSTH_WINDOW, 'YLim', spatial_ylim);
                if i_loc == 1
                    title(ax_raster, event_labels{i_ev}, 'FontSize', 7);
                end
            end

            % X-tick labels (bottom row only)
            if i_loc ~= N_ROWS
                set(ax_psth, 'XTickLabel', []);
            else
                xlabel(ax_psth, 'Time (s)', 'FontSize', 6);
            end

            % Y-label / Y-tick labels (first column only)
            if col == 1
                % define target location label string.
                loc_label = sprintf('Loc %d (%.0f%s)', ...
                    unique_locs(i_loc), ...
                    loc_thetas(i_loc), char(176));

                % make 2-line y-label
                ylabel(ax_psth, {loc_label; 'Mean Spike Rate (sp/s)'}, ...
                    'FontSize', 6);
            else
                set(ax_psth, 'YTickLabel', []);
            end
        end
    end

    %% --- Group 3: Factor Sensitivity (12 cols × 8 rows with rasters) ---
    % Using 8-row grid: odd rows for rasters, even rows for PSTHs
    % Each factor spans 2 rows: row 2k-1 = raster, row 2k = PSTH
    factor_labels = {'Reward', 'Salience', 'Probability', 'Identity'};
    event_labels_full = {'Target On', 'Saccade', 'Reward'};

    % Options for Group 3: 8 rows × 12 columns grid
    opts_g3 = {'LeftMargin', G3_LEFT_MARGIN, 'Width', G3_WIDTH, ...
        'Height', GRID_HEIGHT, 'BottomMargin', GRID_BOTTOM_MARGIN, ...
        'HeightSpacing', 0.008, 'WidthSpacing', G3_WIDTH_SPACING};

    % Track empty panels for legend placement
    legend_ax = [];
    legend_placed = false;

    for i_factor = 1:4
        for i_loc = 1:n_locations
            for i_ev = 1:3
                % Column: (i_loc-1)*3 + i_ev (gives 1-12)
                col = (i_loc - 1) * 3 + i_ev;

                % Calculate grid indices for 8-row × 12-column grid
                % Raster row: 1, 3, 5, 7 for factors 1, 2, 3, 4
                % PSTH row: 2, 4, 6, 8 for factors 1, 2, 3, 4
                raster_row = (i_factor - 1) * 2 + 1;
                psth_row = (i_factor - 1) * 2 + 2;

                % Calculate subplot indices (row-major order)
                raster_idx = (raster_row - 1) * 12 + col;
                psth_idx = (psth_row - 1) * 12 + col;

                % Create raster axes first (top of pair)
                ax_raster = mySubPlot([8, 12, raster_idx], opts_g3{:});

                % Create PSTH axes last so it's current for barStairsFill
                ax_psth = mySubPlot([8, 12, psth_idx], opts_g3{:});

                if ~isempty(factor_psth_data{i_factor, i_loc, i_ev})
                    d = factor_psth_data{i_factor, i_loc, i_ev};

                    % Plot raster with color-coded conditions
                    hold(ax_raster, 'on');

                    % Get trials for each condition
                    cond1_trials = find(d.cond1_mask);
                    cond2_trials = find(d.cond2_mask);

                    % Plot condition 1 trials (High/Face) - blue
                    y_offset = 0;
                    for i_t = 1:length(cond1_trials)
                        trial_idx = cond1_trials(i_t);
                        aligned_spikes = spike_times - ...
                            event_times_all{i_ev}(trial_idx);
                        spikes_in_window = aligned_spikes(...
                            aligned_spikes >= PSTH_WINDOW(1) & ...
                            aligned_spikes <= PSTH_WINDOW(2));
                        if ~isempty(spikes_in_window)
                            plot(ax_raster, spikes_in_window, ...
                                (y_offset + i_t) * ones(size(spikes_in_window)), ...
                                '.', 'Color', COLOR_HIGH, 'MarkerSize', 1);
                        end
                    end

                    % Small gap between condition groups
                    gap_size = 2;
                    y_offset = length(cond1_trials) + gap_size;

                    % Plot condition 2 trials (Low/Non-face) - orange
                    for i_t = 1:length(cond2_trials)
                        trial_idx = cond2_trials(i_t);
                        aligned_spikes = spike_times - ...
                            event_times_all{i_ev}(trial_idx);
                        spikes_in_window = aligned_spikes(...
                            aligned_spikes >= PSTH_WINDOW(1) & ...
                            aligned_spikes <= PSTH_WINDOW(2));
                        if ~isempty(spikes_in_window)
                            plot(ax_raster, spikes_in_window, ...
                                (y_offset + i_t) * ones(size(spikes_in_window)), ...
                                '.', 'Color', COLOR_LOW, 'MarkerSize', 1);
                        end
                    end

                    total_trials = length(cond1_trials) + gap_size + length(cond2_trials);
                    xline(ax_raster, 0, 'k--', 'LineWidth', 0.5);
                    set(ax_raster, 'XLim', PSTH_WINDOW, ...
                        'YLim', [0 total_trials + 1], ...
                        'XTickLabel', [], 'YTick', [], ...
                        'TickDir', 'Out', 'Box', 'Off');

                    % Plot PSTH with two conditions
                    hold(ax_psth, 'on');

                    hBS1 = barStairsFill(d.bin_centers, ...
                        zeros(size(d.rate1)), d.rate1);
                    delete(hBS1(1:2));
                    set(hBS1(3), 'Color', COLOR_HIGH, 'LineWidth', 1);

                    hBS2 = barStairsFill(d.bin_centers, ...
                        zeros(size(d.rate2)), d.rate2);
                    delete(hBS2(1:2));
                    set(hBS2(3), 'Color', COLOR_LOW, 'LineWidth', 1);

                    xline(ax_psth, 0, 'k--', 'LineWidth', 0.5);
                    set(ax_psth, 'XLim', PSTH_WINDOW, 'YLim', factor_ylim, ...
                        'TickDir', 'Out', 'Box', 'Off', 'FontSize', 6);

                    % Column header (top row only)
                    if i_factor == 1
                        title_str = sprintf('L%d:%s', i_loc, ...
                            event_labels_full{i_ev}(1:3));
                        title(ax_raster, title_str, 'FontSize', 6);
                    end
                else
                    set(ax_raster, 'XLim', PSTH_WINDOW, 'YLim', [0 10], ...
                        'XTickLabel', [], 'YTick', [], ...
                        'TickDir', 'Out', 'Box', 'Off');
                    set(ax_psth, 'XLim', PSTH_WINDOW, 'YLim', factor_ylim, ...
                        'TickDir', 'Out', 'Box', 'Off', 'FontSize', 6);

                    % Place legend in Probability row, Location 2 or 4
                    if ~legend_placed && i_factor == 3 && ...
                            (i_loc == 2 || i_loc == 4) && i_ev == 2
                        legend_ax = ax_psth;
                        legend_placed = true;
                    else
                        text(0.5, 0.5, '-', 'Parent', ax_psth, ...
                            'HorizontalAlignment', 'center', ...
                            'FontSize', 6, 'Color', [0.7 0.7 0.7]);
                    end

                    if i_factor == 1
                        title_str = sprintf('L%d:%s', i_loc, ...
                            event_labels_full{i_ev}(1:3));
                        title(ax_raster, title_str, 'FontSize', 6);
                    end
                end

                % X-tick labels (bottom row only)
                if i_factor ~= N_ROWS
                    set(ax_psth, 'XTickLabel', []);
                else
                    xlabel(ax_psth, 'Time (s)', 'FontSize', 5);
                end

                % Y-label / Y-tick labels (first factor column only)
                if col == 1
                    ylabel(ax_psth, {factor_labels{i_factor}; ...
                        'Mean Spike Rate (sp/s)'}, 'FontSize', 6);
                else
                    set(ax_psth, 'YTickLabel', []);
                end
            end
        end
    end

    %% --- Place legend in empty panel ---
    if ~isempty(legend_ax)
        hold(legend_ax, 'on');
        yl = get(legend_ax, 'YLim');
        xl = get(legend_ax, 'XLim');
        y_mid = mean(yl);
        x_start = xl(1) + 0.1 * diff(xl);
        x_end = x_start + 0.3 * diff(xl);
        x_text = x_end + 0.05 * diff(xl);

        plot(legend_ax, [x_start, x_end], [y_mid + 0.2*diff(yl), ...
            y_mid + 0.2*diff(yl)], '-', 'Color', COLOR_HIGH, 'LineWidth', 2);
        text(x_text, y_mid + 0.2*diff(yl), 'High/Face', ...
            'Parent', legend_ax, 'FontSize', 7, 'VerticalAlignment', 'middle');

        plot(legend_ax, [x_start, x_end], [y_mid - 0.2*diff(yl), ...
            y_mid - 0.2*diff(yl)], '-', 'Color', COLOR_LOW, 'LineWidth', 2);
        text(x_text, y_mid - 0.2*diff(yl), 'Low/Non-face', ...
            'Parent', legend_ax, 'FontSize', 7, 'VerticalAlignment', 'middle');
    end

    %% --- Section headers ---
    header_y = GRID_BOTTOM_MARGIN + GRID_HEIGHT + 0.01;

    % Spatial tuning header (above group 2)
    spatial_header_x = G2_LEFT_MARGIN;
    spatial_header_width = G2_WIDTH;
    annotation('textbox', [spatial_header_x, header_y, ...
        spatial_header_width, 0.025], ...
        'String', 'Spatial Tuning', 'FontSize', 8, 'FontWeight', 'bold', ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');

    % Factor sensitivity header (above group 3)
    factor_header_x = G3_LEFT_MARGIN;
    factor_header_width = G3_WIDTH;
    annotation('textbox', [factor_header_x, header_y, ...
        factor_header_width, 0.025], ...
        'String', 'Factor Sensitivity by Location', 'FontSize', 8, ...
        'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center');

    %% --- Final Formatting ---
    all_axes = findall(fig, 'Type', 'Axes');
    set(all_axes, 'TickDir', 'Out', 'LineWidth', 0.5, 'Box', 'Off');

    % if the 'theme' command is one we have, make sure 'theme' is set to
    % 'light':
    if ~isempty(which('theme'))
        theme('light');
    end
    drawnow;

    % Save PDF
    pdfSave(output_filename, figSize / 72, fig);

    fprintf('Saved: %s\n', output_filename);

    close(fig);
end

end
