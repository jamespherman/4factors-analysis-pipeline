function generate_interneuron_psths(varargin)
% GENERATE_INTERNEURON_PSTHS Generate interaction PSTHs for SC interneurons.
%
%   generate_interneuron_psths()
%   generate_interneuron_psths('SessionID', 'Feynman_08_12_2025_SC')
%
%   This function generates two-panel PSTH figures for each SC neuron
%   classified as 'interneuron' by screen_sc_neurons.m.
%
%   Two figures are generated per interneuron:
%     1. Target-aligned PSTHs (aligned to target onset)
%     2. Saccade-aligned PSTHs (aligned to saccade onset)
%
%   Panel 1 (Left): Face vs Non-face, High-salience vs Low-salience
%   Panel 2 (Right): 4 combinations of Reward x Probability
%
%   Trial selection:
%     - Successful gSac_4factors trials only
%     - Trials within the neuron's stable firing rate period
%     - Target at the location with largest |post - pre| activity difference
%
%   OPTIONAL INPUTS (Name-Value pairs):
%       'SessionID'  - Process only this session (default: all SC sessions)
%       'ForceRerun' - Regenerate existing figures (default: false)
%
%   OUTPUT:
%       Figures saved to: figures/interneuron_psths/
%         {session_id}_interneuron_{cluster_id}_targetOn_psth.pdf
%         {session_id}_interneuron_{cluster_id}_saccadeOn_psth.pdf

%% Parse inputs
p = inputParser;
addParameter(p, 'SessionID', '', @ischar);
addParameter(p, 'ForceRerun', false, @islogical);
parse(p, varargin{:});
session_filter = p.Results.SessionID;
force_rerun = p.Results.ForceRerun;

%% Setup paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

%% Create output directory
output_dir = fullfile(project_root, 'figures', 'interneuron_psths');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Load manifest
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
manifest = readtable(manifest_path);

% Filter to SC sessions only
sc_sessions = manifest(strcmp(manifest.brain_area, 'SC'), :);

% If session filter specified, apply it
if ~isempty(session_filter)
    sc_sessions = sc_sessions(strcmp(sc_sessions.unique_id, session_filter), :);
    if isempty(sc_sessions)
        error('Session %s not found or is not an SC session.', session_filter);
    end
end

fprintf('Processing %d SC session(s)...\n', height(sc_sessions));

%% Define colors
% Panel 1 (Face/Salience):
%   Face = teal, Non-face = olive green, High-sal = dark green, Low-sal = gray
% Panel 2 (Reward x Probability):
%   High-rwd/Low-prb = dark purple, High-rwd/High-prb = purple
%   Low-rwd/High-prb = magenta, Low-rwd/Low-prb = gray

colors = struct();
colors.face = [0.0, 0.6, 0.6];      % teal
colors.nonface = [0.5, 0.7, 0.1];   % olive green
colors.high_sal = [0.1, 0.4, 0.1];  % dark green (different from face)
colors.low_sal = [0.6, 0.6, 0.6];   % gray

colors.high_rwd_low_prb = [0.3, 0.0, 0.4];  % dark purple
colors.high_rwd_high_prb = [0.4, 0.2, 0.8]; % purple
colors.low_rwd_high_prb = [0.8, 0.2, 0.6];  % magenta
colors.low_rwd_low_prb = [0.6, 0.6, 0.6];   % gray

%% Analysis parameters
% Target-aligned time window
t_start_target = -0.2;
t_end_target = 0.6;

% Saccade-aligned time window (centered on saccade)
t_start_saccade = -0.4;
t_end_saccade = 0.4;

bin_width = 0.001;  % 1ms bins for smooth traces
smooth_sigma = 0.020;  % 20ms Gaussian smoothing kernel

% Windows for finding best location (relative to target onset)
pre_window = [-0.2, 0];
post_window = [0.05, 0.6];

%% Process each session
codes = initCodes();
total_interneurons = 0;
total_figures = 0;

for i_session = 1:height(sc_sessions)
    session_id = sc_sessions.unique_id{i_session};
    fprintf('\n=== Processing session: %s ===\n', session_id);

    % Load session data
    local_data_path = fullfile(project_root, 'data', 'processed', ...
        session_id, [session_id, '_session_data.mat']);

    if ~exist(local_data_path, 'file')
        fprintf('  Session data not found, skipping.\n');
        continue;
    end

    load(local_data_path, 'session_data');

    % Check for screening info with neuron_class
    if ~isfield(session_data, 'analysis') || ...
       ~isfield(session_data.analysis, 'screening_info')
        fprintf('  No screening_info found, skipping.\n');
        continue;
    end

    screening_info = session_data.analysis.screening_info;

    % screening_info is a struct array - extract neuron_class as cell array
    if ~isfield(screening_info, 'neuron_class')
        fprintf('  No neuron_class field in screening_info, skipping.\n');
        continue;
    end

    neuron_class = {screening_info.neuron_class};
    cluster_ids = [screening_info.cluster_id];

    % Find interneurons
    interneuron_mask = strcmp(neuron_class, 'interneuron');
    interneuron_indices = find(interneuron_mask);
    n_interneurons = length(interneuron_indices);

    fprintf('  Found %d interneurons out of %d total neurons.\n', ...
        n_interneurons, length(neuron_class));

    if n_interneurons == 0
        continue;
    end

    total_interneurons = total_interneurons + n_interneurons;

    % Extract data needed for all interneurons
    trialInfo = session_data.trialInfo;
    eventTimes = session_data.eventTimes;
    all_spike_times = session_data.spikes.times;
    all_spike_clusters = session_data.spikes.clusters;

    % Get all trial indices for gSac_4factors successful trials
    n_trials_total = length(trialInfo.taskCode);
    all_trial_indices = 1:n_trials_total;

    is_4factors = trialInfo.taskCode == codes.uniqueTaskCode_gSac_4factors;
    is_successful = eventTimes.pdsReward > 0;
    base_valid_mask = is_4factors & is_successful;

    if sum(base_valid_mask) == 0
        fprintf('  No valid 4factors trials found, skipping.\n');
        continue;
    end

    % Get event times for ALL trials (we'll filter per-neuron)
    all_target_on_times = eventTimes.targetOn;
    all_saccade_on_times = eventTimes.saccadeOnset;

    % Get condition info for ALL trials
    all_target_locs = trialInfo.targetLocIdx;
    all_reward_dur = trialInfo.rewardDuration;
    all_salience = trialInfo.salience;
    all_stimType = trialInfo.stimType;
    all_block_num = trialInfo.blockNumber;

    %% Process each interneuron
    for i_neuron = 1:n_interneurons
        neuron_idx = interneuron_indices(i_neuron);
        cluster_id = cluster_ids(neuron_idx);

        % Check if BOTH figures already exist
        fig_filename_target = sprintf('%s_interneuron_%03d_targetOn_psth.pdf', ...
            session_id, cluster_id);
        fig_filename_saccade = sprintf('%s_interneuron_%03d_saccadeOn_psth.pdf', ...
            session_id, cluster_id);
        fig_path_target = fullfile(output_dir, fig_filename_target);
        fig_path_saccade = fullfile(output_dir, fig_filename_saccade);

        if exist(fig_path_target, 'file') && exist(fig_path_saccade, 'file') && ~force_rerun
            fprintf('  Neuron %d (cluster %d): Figures exist, skipping.\n', ...
                i_neuron, cluster_id);
            continue;
        end

        fprintf('  Processing neuron %d/%d (cluster %d)...\n', ...
            i_neuron, n_interneurons, cluster_id);

        % Get spike times for this neuron
        spike_times = all_spike_times(all_spike_clusters == cluster_id);

        if isempty(spike_times)
            fprintf('    No spikes found, skipping.\n');
            continue;
        end

        % Get this neuron's stable trial range
        stable_range = screening_info(neuron_idx).stable_trial_range;
        if isempty(stable_range) || ~screening_info(neuron_idx).used_stable_period
            % No stable period detected - use all valid trials
            fprintf('    No stable period, using all valid trials.\n');
            neuron_valid_mask = base_valid_mask;
        else
            % Filter to stable trial range
            stable_start = stable_range(1);
            stable_end = stable_range(2);
            is_in_stable = (all_trial_indices >= stable_start) & (all_trial_indices <= stable_end);
            neuron_valid_mask = base_valid_mask & is_in_stable';
            fprintf('    Stable period: trials %d-%d (%d stable valid trials).\n', ...
                stable_start, stable_end, sum(neuron_valid_mask));
        end

        % Get indices of valid trials for this neuron
        valid_trials = find(neuron_valid_mask);

        if isempty(valid_trials)
            fprintf('    No valid trials in stable period, skipping.\n');
            continue;
        end

        % Extract data for valid trials
        target_locs = all_target_locs(valid_trials);
        unique_locs = unique(target_locs);
        target_on_times = all_target_on_times(valid_trials);
        saccade_on_times = all_saccade_on_times(valid_trials);

        % Condition info for valid trials
        reward_dur = all_reward_dur(valid_trials);
        is_high_reward = reward_dur > 200;
        is_low_reward = reward_dur < 200;

        salience = all_salience(valid_trials);
        is_high_salience = salience == 1;
        is_low_salience = salience == 2;

        stimType = all_stimType(valid_trials);
        is_face = stimType == 1;
        is_nonface = stimType == 2;

        block_num = all_block_num(valid_trials);
        is_high_prob = (block_num == 1 & target_locs == 1) | ...
                       (block_num == 2 & target_locs == 3);
        is_low_prob = (block_num == 1 & ismember(target_locs, [2, 3, 4])) | ...
                      (block_num == 2 & ismember(target_locs, [1, 2, 4]));

        % Find best target location (max absolute pre/post difference)
        % Using target-aligned activity to determine best location
        best_loc = find_best_location(spike_times, target_on_times, ...
            target_locs, unique_locs, pre_window, post_window);

        fprintf('    Best location: %d\n', best_loc);

        % Select trials at best location (same trials for both alignments)
        loc_mask = target_locs == best_loc;
        loc_target_times = target_on_times(loc_mask);
        loc_saccade_times = saccade_on_times(loc_mask);

        % Filter out trials with invalid saccade times (NaN or 0)
        valid_saccade_mask = ~isnan(loc_saccade_times) & loc_saccade_times > 0;

        % Extract condition masks for trials at this location
        loc_is_face = is_face(loc_mask);
        loc_is_nonface = is_nonface(loc_mask);
        loc_is_high_sal = is_high_salience(loc_mask);
        loc_is_low_sal = is_low_salience(loc_mask);
        loc_is_high_rwd = is_high_reward(loc_mask);
        loc_is_low_rwd = is_low_reward(loc_mask);
        loc_is_high_prb = is_high_prob(loc_mask);
        loc_is_low_prb = is_low_prob(loc_mask);

        % Reward x Probability combinations
        hr_lp = loc_is_high_rwd & loc_is_low_prb;
        hr_hp = loc_is_high_rwd & loc_is_high_prb;
        lr_hp = loc_is_low_rwd & loc_is_high_prb;
        lr_lp = loc_is_low_rwd & loc_is_low_prb;

        fprintf('    Trials at location: %d (face=%d, nonface=%d, hi_sal=%d, lo_sal=%d)\n', ...
            sum(loc_mask), sum(loc_is_face), sum(loc_is_nonface), ...
            sum(loc_is_high_sal), sum(loc_is_low_sal));

        % Skip if too few trials at best location
        min_trials_per_location = 10;
        if sum(loc_mask) < min_trials_per_location
            fprintf('    Too few trials at best location (%d < %d), skipping.\n', ...
                sum(loc_mask), min_trials_per_location);
            continue;
        end

        %% Compute firing rates for TARGET-ALIGNED figure
        [time_vec_tgt, fr_face_tgt] = compute_smoothed_fr(spike_times, ...
            loc_target_times(loc_is_face), t_start_target, t_end_target, bin_width, smooth_sigma);
        [~, fr_nonface_tgt] = compute_smoothed_fr(spike_times, ...
            loc_target_times(loc_is_nonface), t_start_target, t_end_target, bin_width, smooth_sigma);
        [~, fr_high_sal_tgt] = compute_smoothed_fr(spike_times, ...
            loc_target_times(loc_is_high_sal), t_start_target, t_end_target, bin_width, smooth_sigma);
        [~, fr_low_sal_tgt] = compute_smoothed_fr(spike_times, ...
            loc_target_times(loc_is_low_sal), t_start_target, t_end_target, bin_width, smooth_sigma);

        [~, fr_hr_lp_tgt] = compute_smoothed_fr(spike_times, ...
            loc_target_times(hr_lp), t_start_target, t_end_target, bin_width, smooth_sigma);
        [~, fr_hr_hp_tgt] = compute_smoothed_fr(spike_times, ...
            loc_target_times(hr_hp), t_start_target, t_end_target, bin_width, smooth_sigma);
        [~, fr_lr_hp_tgt] = compute_smoothed_fr(spike_times, ...
            loc_target_times(lr_hp), t_start_target, t_end_target, bin_width, smooth_sigma);
        [~, fr_lr_lp_tgt] = compute_smoothed_fr(spike_times, ...
            loc_target_times(lr_lp), t_start_target, t_end_target, bin_width, smooth_sigma);

        %% Compute firing rates for SACCADE-ALIGNED figure
        % Use only trials with valid saccade times, but same condition masks
        loc_saccade_valid = loc_saccade_times(valid_saccade_mask);
        face_sac_mask = loc_is_face(valid_saccade_mask);
        nonface_sac_mask = loc_is_nonface(valid_saccade_mask);
        high_sal_sac_mask = loc_is_high_sal(valid_saccade_mask);
        low_sal_sac_mask = loc_is_low_sal(valid_saccade_mask);
        hr_lp_sac_mask = hr_lp(valid_saccade_mask);
        hr_hp_sac_mask = hr_hp(valid_saccade_mask);
        lr_hp_sac_mask = lr_hp(valid_saccade_mask);
        lr_lp_sac_mask = lr_lp(valid_saccade_mask);

        [time_vec_sac, fr_face_sac] = compute_smoothed_fr(spike_times, ...
            loc_saccade_valid(face_sac_mask), t_start_saccade, t_end_saccade, bin_width, smooth_sigma);
        [~, fr_nonface_sac] = compute_smoothed_fr(spike_times, ...
            loc_saccade_valid(nonface_sac_mask), t_start_saccade, t_end_saccade, bin_width, smooth_sigma);
        [~, fr_high_sal_sac] = compute_smoothed_fr(spike_times, ...
            loc_saccade_valid(high_sal_sac_mask), t_start_saccade, t_end_saccade, bin_width, smooth_sigma);
        [~, fr_low_sal_sac] = compute_smoothed_fr(spike_times, ...
            loc_saccade_valid(low_sal_sac_mask), t_start_saccade, t_end_saccade, bin_width, smooth_sigma);

        [~, fr_hr_lp_sac] = compute_smoothed_fr(spike_times, ...
            loc_saccade_valid(hr_lp_sac_mask), t_start_saccade, t_end_saccade, bin_width, smooth_sigma);
        [~, fr_hr_hp_sac] = compute_smoothed_fr(spike_times, ...
            loc_saccade_valid(hr_hp_sac_mask), t_start_saccade, t_end_saccade, bin_width, smooth_sigma);
        [~, fr_lr_hp_sac] = compute_smoothed_fr(spike_times, ...
            loc_saccade_valid(lr_hp_sac_mask), t_start_saccade, t_end_saccade, bin_width, smooth_sigma);
        [~, fr_lr_lp_sac] = compute_smoothed_fr(spike_times, ...
            loc_saccade_valid(lr_lp_sac_mask), t_start_saccade, t_end_saccade, bin_width, smooth_sigma);

        %% Determine shared y-axis limits (same for both figures)
        all_frs = [fr_face_tgt(:); fr_nonface_tgt(:); fr_high_sal_tgt(:); fr_low_sal_tgt(:); ...
                   fr_hr_lp_tgt(:); fr_hr_hp_tgt(:); fr_lr_hp_tgt(:); fr_lr_lp_tgt(:); ...
                   fr_face_sac(:); fr_nonface_sac(:); fr_high_sal_sac(:); fr_low_sal_sac(:); ...
                   fr_hr_lp_sac(:); fr_hr_hp_sac(:); fr_lr_hp_sac(:); fr_lr_lp_sac(:)];
        all_frs = all_frs(~isnan(all_frs));

        if isempty(all_frs) || max(all_frs) == 0
            fprintf('    No valid firing rate data, skipping.\n');
            continue;
        end

        y_max_raw = max(all_frs);
        % Round up to nice number for y-axis
        y_max = ceil(y_max_raw / 10) * 10;
        if y_max < 10
            y_max = ceil(y_max_raw * 2) / 2;  % Round to nearest 0.5
        end
        y_max = max(y_max, 5);  % Minimum y_max of 5

        %% Create and save TARGET-ALIGNED figure
        plot_and_save_psth_figure(fig_path_target, time_vec_tgt, ...
            fr_face_tgt, fr_nonface_tgt, fr_high_sal_tgt, fr_low_sal_tgt, ...
            fr_hr_lp_tgt, fr_hr_hp_tgt, fr_lr_hp_tgt, fr_lr_lp_tgt, ...
            loc_is_face, loc_is_nonface, loc_is_high_sal, loc_is_low_sal, ...
            hr_lp, hr_hp, lr_hp, lr_lp, ...
            colors, t_start_target, t_end_target, y_max);
        fprintf('    Saved: %s\n', fig_filename_target);
        total_figures = total_figures + 1;

        %% Create and save SACCADE-ALIGNED figure
        plot_and_save_psth_figure(fig_path_saccade, time_vec_sac, ...
            fr_face_sac, fr_nonface_sac, fr_high_sal_sac, fr_low_sal_sac, ...
            fr_hr_lp_sac, fr_hr_hp_sac, fr_lr_hp_sac, fr_lr_lp_sac, ...
            face_sac_mask, nonface_sac_mask, high_sal_sac_mask, low_sal_sac_mask, ...
            hr_lp_sac_mask, hr_hp_sac_mask, lr_hp_sac_mask, lr_lp_sac_mask, ...
            colors, t_start_saccade, t_end_saccade, y_max);
        fprintf('    Saved: %s\n', fig_filename_saccade);
        total_figures = total_figures + 1;
    end
end

fprintf('\n=== Complete ===\n');
fprintf('Total interneurons processed: %d\n', total_interneurons);
fprintf('Total figures generated: %d\n', total_figures);
fprintf('Figures saved to: %s\n', output_dir);

end

%% Helper Functions

function plot_and_save_psth_figure(fig_path, time_vec, ...
    fr_face, fr_nonface, fr_high_sal, fr_low_sal, ...
    fr_hr_lp, fr_hr_hp, fr_lr_hp, fr_lr_lp, ...
    mask_face, mask_nonface, mask_high_sal, mask_low_sal, ...
    mask_hr_lp, mask_hr_hp, mask_lr_hp, mask_lr_lp, ...
    colors, t_start, t_end, y_max)
% Create and save a two-panel PSTH figure.

fig = figure('Color', 'w', 'Position', [100, 100, 192, 96], ...
    'MenuBar', 'none', 'ToolBar', 'none');

% Panel 1: Face/Salience
ax1 = subplot(1, 2, 1);
hold on;

% Plot traces (order: face, nonface, high_sal, low_sal)
if any(mask_face)
    plot(time_vec, fr_face, 'Color', colors.face, 'LineWidth', 0.75);
end
if any(mask_nonface)
    plot(time_vec, fr_nonface, 'Color', colors.nonface, 'LineWidth', 0.75);
end
if any(mask_high_sal)
    plot(time_vec, fr_high_sal, 'Color', colors.high_sal, 'LineWidth', 0.75);
end
if any(mask_low_sal)
    plot(time_vec, fr_low_sal, 'Color', colors.low_sal, 'LineWidth', 0.75);
end

xline(0, '--', 'Color', [0.5, 0.5, 0.5]);
xlim([t_start, t_end]);
ylim([0, y_max]);

% Add y-tick labels before drawAxes removes native axes
% Position labels relative to plot start
y_mid = round(y_max / 2);
y_ticks = [0, y_mid, y_max];
label_x = t_start - 0.15;  % Position labels to the left of plot
for yt = y_ticks
    text(label_x, yt, num2str(yt), 'FontSize', 5, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

% Draw disconnected axes (matching user's example format)
% xAxisPos is slightly to the left of t_start to create gap
x_params = [t_start - 0.05, t_start, t_end];
y_params = [-y_max * 0.07, 0, y_max];
drawAxes(ax1, x_params, y_params, 0.06, 5, 3);

% Panel 2: Reward x Probability
ax2 = subplot(1, 2, 2);
hold on;

% Plot traces
if any(mask_hr_lp)
    plot(time_vec, fr_hr_lp, 'Color', colors.high_rwd_low_prb, 'LineWidth', 0.75);
end
if any(mask_hr_hp)
    plot(time_vec, fr_hr_hp, 'Color', colors.high_rwd_high_prb, 'LineWidth', 0.75);
end
if any(mask_lr_hp)
    plot(time_vec, fr_lr_hp, 'Color', colors.low_rwd_high_prb, 'LineWidth', 0.75);
end
if any(mask_lr_lp)
    plot(time_vec, fr_lr_lp, 'Color', colors.low_rwd_low_prb, 'LineWidth', 0.75);
end

xline(0, '--', 'Color', [0.5, 0.5, 0.5]);
xlim([t_start, t_end]);
ylim([0, y_max]);

% Draw disconnected axes (same params as ax1)
drawAxes(ax2, x_params, y_params, 0.06, 5, 3);

% Apply light theme if available
try
    theme('light');
catch
end

% Set renderer to Painters for true vector graphics
set(fig, 'Renderer', 'painters');

% Save figure as vector PDF
set(fig, 'PaperUnits', 'Inches', 'PaperSize', fig.Position(3:4)/72);
set(fig, 'PaperPositionMode', 'auto');
print(fig, fig_path, '-dpdf', '-painters');
close(fig);

end

function best_loc = find_best_location(spike_times, event_times, ...
    trial_locs, unique_locs, pre_window, post_window)
% Find target location with largest absolute pre/post activity difference.

max_diff = -inf;
best_loc = unique_locs(1);

for loc = unique_locs(:)'
    loc_mask = trial_locs == loc;
    loc_events = event_times(loc_mask);

    if isempty(loc_events)
        continue;
    end

    % Compute mean firing rate in pre and post windows
    pre_spikes = 0;
    post_spikes = 0;
    pre_duration = diff(pre_window);
    post_duration = diff(post_window);

    for i_trial = 1:length(loc_events)
        t0 = loc_events(i_trial);
        pre_spikes = pre_spikes + sum(spike_times >= t0 + pre_window(1) & ...
                                       spike_times < t0 + pre_window(2));
        post_spikes = post_spikes + sum(spike_times >= t0 + post_window(1) & ...
                                         spike_times < t0 + post_window(2));
    end

    n_trials = length(loc_events);
    pre_fr = pre_spikes / (n_trials * pre_duration);
    post_fr = post_spikes / (n_trials * post_duration);

    abs_diff = abs(post_fr - pre_fr);

    if abs_diff > max_diff
        max_diff = abs_diff;
        best_loc = loc;
    end
end

end

function [time_vec, smoothed_fr] = compute_smoothed_fr(spike_times, ...
    event_times, t_start, t_end, bin_width, smooth_sigma)
% Compute smoothed firing rate aligned to events.

% Default time vector
time_vec = t_start:bin_width:t_end;

if isempty(event_times)
    smoothed_fr = nan(1, length(time_vec));
    return;
end

% Bin spikes using alignAndBinSpikes
[binned_counts, bin_centers, ~] = alignAndBinSpikes(spike_times, ...
    event_times, t_start, t_end, bin_width);

% Average across trials
mean_counts = mean(binned_counts, 1, 'omitnan');

% Ensure row vector
mean_counts = mean_counts(:)';

% Convert to firing rate (spikes/second)
mean_fr = mean_counts / bin_width;

% Gaussian smoothing
if smooth_sigma > 0
    % Create Gaussian kernel
    kernel_width = ceil(4 * smooth_sigma / bin_width);
    kernel_x = -kernel_width:kernel_width;
    kernel = exp(-0.5 * (kernel_x * bin_width / smooth_sigma).^2);
    kernel = kernel / sum(kernel);

    % Smooth the firing rate
    smoothed_fr = conv(mean_fr, kernel, 'same');
else
    smoothed_fr = mean_fr;
end

% Ensure output is row vector
smoothed_fr = smoothed_fr(:)';
time_vec = bin_centers(:)';

end
