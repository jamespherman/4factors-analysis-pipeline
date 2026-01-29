%% plot_behavioral_summary_pooled.m
%
% Generates a summary figure of behavioral metrics pooled across all
% subjects and sessions. Displays median reaction time, endpoint error,
% and saccadic peak velocity for correct trials, broken down by 16
% conditions (8 factor levels x 2 target locations).
%
% The 8 factor levels are:
%   1. High Salience (bullseye)
%   2. Low Salience (bullseye)
%   3. Face
%   4. Non-Face
%   5. High Reward x Low Probability
%   6. High Reward x High Probability
%   7. Low Reward x High Probability
%   8. Low Reward x Low Probability
%
% The 2 target locations are: Location 1 and Location 3
%
% Author: Claude
% Date: 2026-01-29

%% Setup
clear; clc; close all;

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

%% Configuration
% Bootstrap parameters for 95% CI
n_bootstrap = 1000;
ci_alpha = 0.05;

% Target locations to include
target_locations = [1, 3];

% Condition labels and definitions
% Each condition specifies: {name, filter_function}
% The filter function takes trialInfo and returns a logical mask
condition_defs = {
    'High Sal',     @(ti, et) filter_high_salience(ti);
    'Low Sal',      @(ti, et) filter_low_salience(ti);
    'Face',         @(ti, et) filter_face(ti);
    'Non-Face',     @(ti, et) filter_nonface(ti);
    'HiRwd-LoProb', @(ti, et) filter_hirwd_loprob(ti);
    'HiRwd-HiProb', @(ti, et) filter_hirwd_hiprob(ti);
    'LoRwd-HiProb', @(ti, et) filter_lorwd_hiprob(ti);
    'LoRwd-LoProb', @(ti, et) filter_lorwd_loprob(ti);
};

n_conditions = size(condition_defs, 1);
n_locations = length(target_locations);
n_total_bars = n_conditions * n_locations;

% Main effect condition definitions (for Figure 3)
% These are the 4 factors x 2 levels = 8 main effects (not combinations)
main_effect_defs = {
    'High Sal',   @(ti, et) filter_high_salience(ti);
    'Low Sal',    @(ti, et) filter_low_salience(ti);
    'Face',       @(ti, et) filter_face(ti);
    'Non-Face',   @(ti, et) filter_nonface(ti);
    'High Rwd',   @(ti, et) filter_high_reward(ti);
    'Low Rwd',    @(ti, et) filter_low_reward(ti);
    'High Prob',  @(ti, et) filter_high_probability(ti);
    'Low Prob',   @(ti, et) filter_low_probability(ti);
};
n_main_effects = size(main_effect_defs, 1);

%% Load Manifest
manifest_path = fullfile(project_root, 'config', 'session_manifest.csv');
if ~exist(manifest_path, 'file')
    error('Session manifest not found at: %s', manifest_path);
end
manifest = readtable(manifest_path);
fprintf('Found %d sessions in manifest.\n', height(manifest));

%% Initialize Data Storage
% Store all trial-level data for pooling (Figure 1 & 2 conditions)
all_RT = [];
all_EPE = [];  % Endpoint Error
all_PV = [];   % Peak Velocity
all_condition_idx = [];
all_location_idx = [];
all_session_idx = [];

% Store data for main effects (Figure 3 conditions)
all_RT_main = [];
all_EPE_main = [];
all_PV_main = [];
all_main_effect_idx = [];
all_location_idx_main = [];
all_session_idx_main = [];

%% Load and Process Each Session
codes = initCodes();

for i_session = 1:height(manifest)
    session_id = manifest.unique_id{i_session};
    fprintf('Processing session %d/%d: %s\n', i_session, height(manifest), session_id);

    % Try local processed directory first
    local_processed_dir = fullfile(project_root, 'data', 'processed', session_id);
    local_session_data_path = fullfile(local_processed_dir, [session_id, '_session_data.mat']);

    if exist(local_session_data_path, 'file')
        data_path = local_session_data_path;
    else
        % Fall back to OneDrive
        one_drive_path = findOneDrive();
        data_path = fullfile(one_drive_path, 'Neuronal Data Analysis', ...
            session_id, [session_id '_session_data.mat']);
    end

    if ~exist(data_path, 'file')
        warning('Data not found for %s. Skipping.', session_id);
        continue;
    end

    load(data_path, 'session_data');

    trialInfo = session_data.trialInfo;
    eventTimes = session_data.eventTimes;

    % Determine number of trials (trialInfo can be a table or struct)
    if istable(trialInfo)
        n_trials_total = height(trialInfo);
    else
        % It's a struct - get length from any field
        n_trials_total = length(trialInfo.taskCode);
    end

    % --- Identify gSac_4factors trials ---
    is_gsac_4factors = trialInfo.taskCode == codes.uniqueTaskCode_gSac_4factors;

    % --- Identify correct trials (reward delivered) ---
    % Check for non-empty reward times in eventTimes
    if isfield(eventTimes, 'rewardCell')
        is_correct = ~cellfun(@isempty, eventTimes.rewardCell);
    elseif isfield(eventTimes, 'reward')
        is_correct = eventTimes.reward > 0;
    else
        % Fall back to checking pdsReward
        is_correct = eventTimes.pdsReward > 0;
    end

    % Master mask: correct gSac_4factors trials
    master_mask = is_gsac_4factors & is_correct;
    n_valid_trials = sum(master_mask);

    if n_valid_trials == 0
        fprintf('  No valid trials found. Skipping.\n');
        continue;
    end

    fprintf('  Found %d correct gSac_4factors trials.\n', n_valid_trials);

    % --- Extract behavioral metrics ---
    % Reaction Time: compute from event times (SRT field may contain timestamps)
    % RT = time from fixOff (go signal) to saccadeOnset
    sac_onset = eventTimes.pdsSaccadeOnset(:);
    fix_off = eventTimes.pdsFixOff(:);

    % Replace invalid values (-1 or NaN) with NaN
    sac_onset(sac_onset < 0) = NaN;
    fix_off(fix_off < 0) = NaN;

    RT = sac_onset - fix_off;

    % Filter out unreasonable RT values (< 0 or > 1 second)
    RT(RT < 0 | RT > 1) = NaN;

    % Endpoint Error
    if isfield(trialInfo, 'postSacXY')
        post_sac_xy = trialInfo.postSacXY;
        targ_x = trialInfo.targDegX;
        targ_y = trialInfo.targDegY;

        % Compute Euclidean distance
        if size(post_sac_xy, 2) >= 2
            dx = post_sac_xy(:, 1) - targ_x(:);
            dy = post_sac_xy(:, 2) - targ_y(:);
            EPE = sqrt(dx.^2 + dy.^2);
        else
            EPE = nan(n_trials_total, 1);
        end
    else
        EPE = nan(n_trials_total, 1);
    end

    % Peak Velocity
    if isfield(trialInfo, 'peakVel')
        PV = trialInfo.peakVel;
    else
        PV = nan(n_trials_total, 1);
    end
    PV = PV(:);

    % --- Extract location indices ---
    target_loc_idx = trialInfo.targetLocIdx;

    % --- Classify trials into conditions ---
    % Each condition is computed independently (trials can match multiple conditions)
    % We store boolean masks per condition for each trial
    condition_masks = false(n_trials_total, n_conditions);

    for i_cond = 1:n_conditions
        filter_fn = condition_defs{i_cond, 2};
        cond_mask = filter_fn(trialInfo, eventTimes);
        cond_mask = logical(cond_mask(:));

        if length(cond_mask) ~= n_trials_total
            warning('Condition mask size mismatch for condition %d: got %d, expected %d', ...
                i_cond, length(cond_mask), n_trials_total);
            continue;
        end

        condition_masks(:, i_cond) = cond_mask;
    end

    % Assign location index (1 or 2 for locations 1 and 3)
    trial_location = zeros(n_trials_total, 1);
    for i_loc = 1:n_locations
        loc_mask = (target_loc_idx == target_locations(i_loc));
        loc_mask = logical(loc_mask(:));
        trial_location(loc_mask) = i_loc;
    end

    % --- Append data for each condition x location combination ---
    % Trials can contribute to multiple conditions (since conditions are independent)
    for i_cond = 1:n_conditions
        for i_loc = 1:n_locations
            % Valid trials: correct gSac_4factors, matches this condition, matches this location
            valid_mask = master_mask & condition_masks(:, i_cond) & (trial_location == i_loc);
            n_valid = sum(valid_mask);

            if n_valid > 0
                all_RT = [all_RT; RT(valid_mask)];
                all_EPE = [all_EPE; EPE(valid_mask)];
                all_PV = [all_PV; PV(valid_mask)];
                all_condition_idx = [all_condition_idx; repmat(i_cond, n_valid, 1)];
                all_location_idx = [all_location_idx; repmat(i_loc, n_valid, 1)];
                all_session_idx = [all_session_idx; repmat(i_session, n_valid, 1)];
            end
        end
    end

    % --- Collect data for main effect conditions (Figure 3) ---
    main_effect_masks = false(n_trials_total, n_main_effects);
    for i_me = 1:n_main_effects
        filter_fn = main_effect_defs{i_me, 2};
        me_mask = filter_fn(trialInfo, eventTimes);
        me_mask = logical(me_mask(:));
        if length(me_mask) == n_trials_total
            main_effect_masks(:, i_me) = me_mask;
        end
    end

    % Append main effect data
    for i_me = 1:n_main_effects
        for i_loc = 1:n_locations
            valid_mask = master_mask & main_effect_masks(:, i_me) & (trial_location == i_loc);
            n_valid = sum(valid_mask);

            if n_valid > 0
                all_RT_main = [all_RT_main; RT(valid_mask)];
                all_EPE_main = [all_EPE_main; EPE(valid_mask)];
                all_PV_main = [all_PV_main; PV(valid_mask)];
                all_main_effect_idx = [all_main_effect_idx; repmat(i_me, n_valid, 1)];
                all_location_idx_main = [all_location_idx_main; repmat(i_loc, n_valid, 1)];
                all_session_idx_main = [all_session_idx_main; repmat(i_session, n_valid, 1)];
            end
        end
    end
end

fprintf('\nTotal pooled trials: %d\n', length(all_RT));

%% Compute Medians and 95% CIs for Each Condition x Location
medians_RT = nan(n_conditions, n_locations);
ci_RT = nan(n_conditions, n_locations, 2);

medians_EPE = nan(n_conditions, n_locations);
ci_EPE = nan(n_conditions, n_locations, 2);

medians_PV = nan(n_conditions, n_locations);
ci_PV = nan(n_conditions, n_locations, 2);

for i_cond = 1:n_conditions
    for i_loc = 1:n_locations
        mask = (all_condition_idx == i_cond) & (all_location_idx == i_loc);
        n_trials = sum(mask);

        fprintf('Condition %d (%s), Location %d: %d trials\n', ...
            i_cond, condition_defs{i_cond, 1}, target_locations(i_loc), n_trials);

        if n_trials < 5
            continue;
        end

        % Reaction Time
        rt_data = all_RT(mask);
        rt_data = rt_data(~isnan(rt_data));
        if ~isempty(rt_data)
            medians_RT(i_cond, i_loc) = median(rt_data);
            ci_RT(i_cond, i_loc, :) = bootstrap_ci(rt_data, n_bootstrap, ci_alpha);
        end

        % Endpoint Error
        epe_data = all_EPE(mask);
        epe_data = epe_data(~isnan(epe_data));
        if ~isempty(epe_data)
            medians_EPE(i_cond, i_loc) = median(epe_data);
            ci_EPE(i_cond, i_loc, :) = bootstrap_ci(epe_data, n_bootstrap, ci_alpha);
        end

        % Peak Velocity
        pv_data = all_PV(mask);
        pv_data = pv_data(~isnan(pv_data));
        if ~isempty(pv_data)
            medians_PV(i_cond, i_loc) = median(pv_data);
            ci_PV(i_cond, i_loc, :) = bootstrap_ci(pv_data, n_bootstrap, ci_alpha);
        end
    end
end

%% Create Figure
fig = figure('Position', [100, 100, 600, 900], 'Color', 'w');

% Colors for the two locations
colors = [0.2 0.4 0.8;   % Location 1 - Blue
          0.8 0.4 0.2];  % Location 3 - Orange

% Bar width
bar_width = 0.15;

% X positions for bars
x_base = 1:n_conditions;

% Metrics to plot
metrics = struct();
metrics(1).name = 'Median Reaction Time (s)';
metrics(1).data = medians_RT;
metrics(1).ci = ci_RT;

metrics(2).name = 'Endpoint Error (deg)';
metrics(2).data = medians_EPE;
metrics(2).ci = ci_EPE;

metrics(3).name = 'Peak Velocity (deg/s)';
metrics(3).data = medians_PV;
metrics(3).ci = ci_PV;

%% Plot Each Row
for i_row = 1:3
    ax = mySubPlot([3, 1, i_row], 'LeftMargin', 0.125);
    hold on;

    y_data = metrics(i_row).data;
    y_ci = metrics(i_row).ci;

    % Plot bars for each location
    for i_loc = 1:n_locations
        x_offset = (i_loc - 1.5) * bar_width * 2;
        x_pos = x_base + x_offset;

        y_vals = y_data(:, i_loc);
        y_lo = squeeze(y_ci(:, i_loc, 1));
        y_hi = squeeze(y_ci(:, i_loc, 2));

        % Use mybarerr for each bar
        for i_cond = 1:n_conditions
            if ~isnan(y_vals(i_cond))
                mybarerr(x_pos(i_cond), y_vals(i_cond), ...
                    [y_lo(i_cond), y_hi(i_cond)], [], colors(i_loc, :), bar_width*0.7);
            end
        end
    end

    % Formatting
    ylabel(metrics(i_row).name, 'FontSize', 11);
    set(ax, 'XTick', x_base, 'XTickLabel', condition_defs(:, 1), ...
        'FontSize', 9, 'TickDir', 'Out');
    xlim([0.5, n_conditions + 0.5]);
    if i_row < 3
        set(ax, 'XColor', 'none')
    end

    % Legend for first row
    if i_row == 1
        legend_entries = arrayfun(@(x) sprintf('Location %d', x), target_locations, ...
            'UniformOutput', false);
        % Create dummy patches for legend
        h_leg = [];
        for i_loc = 1:n_locations
            h_leg(i_loc) = patch(nan, nan, colors(i_loc, :), 'EdgeColor', 'none');
        end
        legend(h_leg, legend_entries, 'Location', 'northeast', 'FontSize', 9);
    end

    hold off;
end

%% Save Figure
figures_dir = fullfile(project_root, 'figures', 'behavior');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

fig_filename = fullfile(figures_dir, 'pooled_behavioral_summary.pdf');
pdfSave(fig_filename, fig.Position(3:4)/72, fig);
fprintf('\nFigure saved to: %s\n', fig_filename);

%% ========================================================================
%  FIGURE 2: Collapsed Across Locations
%  ========================================================================
%  Uses stratified bootstrap to compute CIs that reflect within-location
%  variability rather than across-location variability. For each bootstrap
%  iteration, we resample within each location, compute location-specific
%  medians, then average them. This treats location as a blocking factor.

fprintf('\n--- Computing collapsed statistics (stratified bootstrap) ---\n');

medians_RT_collapsed = nan(n_conditions, 1);
ci_RT_collapsed = nan(n_conditions, 2);

medians_EPE_collapsed = nan(n_conditions, 1);
ci_EPE_collapsed = nan(n_conditions, 2);

medians_PV_collapsed = nan(n_conditions, 1);
ci_PV_collapsed = nan(n_conditions, 2);

for i_cond = 1:n_conditions
    % Get data for each location
    mask_loc1 = (all_condition_idx == i_cond) & (all_location_idx == 1);
    mask_loc2 = (all_condition_idx == i_cond) & (all_location_idx == 2);

    % --- Reaction Time ---
    rt_loc1 = all_RT(mask_loc1);
    rt_loc1 = rt_loc1(~isnan(rt_loc1));
    rt_loc2 = all_RT(mask_loc2);
    rt_loc2 = rt_loc2(~isnan(rt_loc2));

    if ~isempty(rt_loc1) && ~isempty(rt_loc2)
        % Point estimate: average of location-specific medians
        medians_RT_collapsed(i_cond) = (median(rt_loc1) + median(rt_loc2)) / 2;
        % Stratified bootstrap CI
        ci_RT_collapsed(i_cond, :) = stratified_bootstrap_ci(...
            rt_loc1, rt_loc2, n_bootstrap, ci_alpha);
    end

    % --- Endpoint Error ---
    epe_loc1 = all_EPE(mask_loc1);
    epe_loc1 = epe_loc1(~isnan(epe_loc1));
    epe_loc2 = all_EPE(mask_loc2);
    epe_loc2 = epe_loc2(~isnan(epe_loc2));

    if ~isempty(epe_loc1) && ~isempty(epe_loc2)
        medians_EPE_collapsed(i_cond) = (median(epe_loc1) + median(epe_loc2)) / 2;
        ci_EPE_collapsed(i_cond, :) = stratified_bootstrap_ci(...
            epe_loc1, epe_loc2, n_bootstrap, ci_alpha);
    end

    % --- Peak Velocity ---
    pv_loc1 = all_PV(mask_loc1);
    pv_loc1 = pv_loc1(~isnan(pv_loc1));
    pv_loc2 = all_PV(mask_loc2);
    pv_loc2 = pv_loc2(~isnan(pv_loc2));

    if ~isempty(pv_loc1) && ~isempty(pv_loc2)
        medians_PV_collapsed(i_cond) = (median(pv_loc1) + median(pv_loc2)) / 2;
        ci_PV_collapsed(i_cond, :) = stratified_bootstrap_ci(...
            pv_loc1, pv_loc2, n_bootstrap, ci_alpha);
    end

    fprintf('Condition %d (%s) collapsed: RT=%.3f, EPE=%.2f, PV=%.1f\n', ...
        i_cond, condition_defs{i_cond, 1}, ...
        medians_RT_collapsed(i_cond), ...
        medians_EPE_collapsed(i_cond), ...
        medians_PV_collapsed(i_cond));
end

%% Create Collapsed Figure (Figure 2)
fig2 = figure('Position', [100, 100, 600, 900], 'Color', 'w');

% Alternating colors for paired conditions (green for "high/first", purple for "low/second")
color_green = [0.2 0.6 0.3];   % High/first levels
color_purple = [0.5 0.2 0.6];  % Low/second levels
% Color assignment: odd indices (1,3,5,7) -> green, even indices (2,4,6,8) -> purple
collapsed_colors = repmat([color_green; color_purple], 4, 1);

% Bar width for collapsed figure
bar_width_collapsed = 0.25;

% Metrics to plot (collapsed)
metrics_collapsed = struct();
metrics_collapsed(1).name = 'Median Reaction Time (s)';
metrics_collapsed(1).data = medians_RT_collapsed;
metrics_collapsed(1).ci = ci_RT_collapsed;

metrics_collapsed(2).name = 'Endpoint Error (deg)';
metrics_collapsed(2).data = medians_EPE_collapsed;
metrics_collapsed(2).ci = ci_EPE_collapsed;

metrics_collapsed(3).name = 'Peak Velocity (deg/s)';
metrics_collapsed(3).data = medians_PV_collapsed;
metrics_collapsed(3).ci = ci_PV_collapsed;

%% Plot Each Row (Collapsed - Figure 2)
for i_row = 1:3
    ax = mySubPlot([3, 1, i_row], 'LeftMargin', 0.125);
    hold on;

    y_data = metrics_collapsed(i_row).data;
    y_ci = metrics_collapsed(i_row).ci;

    y_lo = y_ci(:, 1);
    y_hi = y_ci(:, 2);

    % Plot bars for each condition with alternating colors
    for i_cond = 1:n_conditions
        if ~isnan(y_data(i_cond))
            mybarerr(i_cond, y_data(i_cond), ...
                [y_lo(i_cond), y_hi(i_cond)], [], collapsed_colors(i_cond, :), bar_width_collapsed);
        end
    end

    % Formatting
    ylabel(metrics_collapsed(i_row).name, 'FontSize', 11);
    set(ax, 'XTick', x_base, 'XTickLabel', condition_defs(:, 1), ...
        'FontSize', 9, 'TickDir', 'Out');
    xlim([0.5, n_conditions + 0.5]);
    if i_row < 3
        set(ax, 'XColor', 'none')
    end

    hold off;
end

%% Save Collapsed Figure
fig2_filename = fullfile(figures_dir, 'pooled_behavioral_summary_collapsed.pdf');
pdfSave(fig2_filename, fig2.Position(3:4)/72, fig2);
fprintf('\nCollapsed figure saved to: %s\n', fig2_filename);

%% ========================================================================
%  FIGURE 3: Main Effects (Collapsed Across Locations)
%  ========================================================================
%  Shows the 8 main factor levels (not combinations):
%  High/Low Salience, Face/Non-Face, High/Low Reward, High/Low Probability

fprintf('\n--- Computing main effect statistics (stratified bootstrap) ---\n');

medians_RT_main = nan(n_main_effects, 1);
ci_RT_main = nan(n_main_effects, 2);

medians_EPE_main = nan(n_main_effects, 1);
ci_EPE_main = nan(n_main_effects, 2);

medians_PV_main = nan(n_main_effects, 1);
ci_PV_main = nan(n_main_effects, 2);

for i_me = 1:n_main_effects
    % Get data for each location
    mask_loc1 = (all_main_effect_idx == i_me) & (all_location_idx_main == 1);
    mask_loc2 = (all_main_effect_idx == i_me) & (all_location_idx_main == 2);

    % --- Reaction Time ---
    rt_loc1 = all_RT_main(mask_loc1);
    rt_loc1 = rt_loc1(~isnan(rt_loc1));
    rt_loc2 = all_RT_main(mask_loc2);
    rt_loc2 = rt_loc2(~isnan(rt_loc2));

    if ~isempty(rt_loc1) && ~isempty(rt_loc2)
        medians_RT_main(i_me) = (median(rt_loc1) + median(rt_loc2)) / 2;
        ci_RT_main(i_me, :) = stratified_bootstrap_ci(rt_loc1, rt_loc2, n_bootstrap, ci_alpha);
    end

    % --- Endpoint Error ---
    epe_loc1 = all_EPE_main(mask_loc1);
    epe_loc1 = epe_loc1(~isnan(epe_loc1));
    epe_loc2 = all_EPE_main(mask_loc2);
    epe_loc2 = epe_loc2(~isnan(epe_loc2));

    if ~isempty(epe_loc1) && ~isempty(epe_loc2)
        medians_EPE_main(i_me) = (median(epe_loc1) + median(epe_loc2)) / 2;
        ci_EPE_main(i_me, :) = stratified_bootstrap_ci(epe_loc1, epe_loc2, n_bootstrap, ci_alpha);
    end

    % --- Peak Velocity ---
    pv_loc1 = all_PV_main(mask_loc1);
    pv_loc1 = pv_loc1(~isnan(pv_loc1));
    pv_loc2 = all_PV_main(mask_loc2);
    pv_loc2 = pv_loc2(~isnan(pv_loc2));

    if ~isempty(pv_loc1) && ~isempty(pv_loc2)
        medians_PV_main(i_me) = (median(pv_loc1) + median(pv_loc2)) / 2;
        ci_PV_main(i_me, :) = stratified_bootstrap_ci(pv_loc1, pv_loc2, n_bootstrap, ci_alpha);
    end

    fprintf('Main effect %d (%s): RT=%.3f, EPE=%.2f, PV=%.1f\n', ...
        i_me, main_effect_defs{i_me, 1}, ...
        medians_RT_main(i_me), medians_EPE_main(i_me), medians_PV_main(i_me));
end

%% Create Main Effects Figure (Figure 3)
fig3 = figure('Position', [100, 100, 600, 900], 'Color', 'w');

% Alternating colors for paired conditions (green for "high/first", purple for "low/second")
main_effect_colors = repmat([color_green; color_purple], 4, 1);

% Bar width
bar_width_main = 0.25;

% Metrics to plot
metrics_main = struct();
metrics_main(1).name = 'Median Reaction Time (s)';
metrics_main(1).data = medians_RT_main;
metrics_main(1).ci = ci_RT_main;

metrics_main(2).name = 'Endpoint Error (deg)';
metrics_main(2).data = medians_EPE_main;
metrics_main(2).ci = ci_EPE_main;

metrics_main(3).name = 'Peak Velocity (deg/s)';
metrics_main(3).data = medians_PV_main;
metrics_main(3).ci = ci_PV_main;

x_base_main = 1:n_main_effects;

%% Plot Each Row (Main Effects - Figure 3)
for i_row = 1:3
    ax = mySubPlot([3, 1, i_row], 'LeftMargin', 0.125);
    hold on;

    y_data = metrics_main(i_row).data;
    y_ci = metrics_main(i_row).ci;

    y_lo = y_ci(:, 1);
    y_hi = y_ci(:, 2);

    % Plot bars for each main effect with alternating colors
    for i_me = 1:n_main_effects
        if ~isnan(y_data(i_me))
            mybarerr(i_me, y_data(i_me), ...
                [y_lo(i_me), y_hi(i_me)], [], main_effect_colors(i_me, :), bar_width_main);
        end
    end

    % Formatting
    ylabel(metrics_main(i_row).name, 'FontSize', 11);
    set(ax, 'XTick', x_base_main, 'XTickLabel', main_effect_defs(:, 1), ...
        'FontSize', 9, 'TickDir', 'Out');
    xlim([0.5, n_main_effects + 0.5]);
    if i_row < 3
        set(ax, 'XColor', 'none')
    end

    hold off;
end

%% Save Main Effects Figure
fig3_filename = fullfile(figures_dir, 'pooled_behavioral_summary_main_effects.pdf');
pdfSave(fig3_filename, fig3.Position(3:4)/72, fig3);
fprintf('\nMain effects figure saved to: %s\n', fig3_filename);

%% Helper Functions

function ci = bootstrap_ci(data, n_boot, alpha)
    % Compute bootstrap confidence interval for the median
    n = length(data);
    boot_medians = nan(n_boot, 1);

    for i = 1:n_boot
        boot_sample = data(randi(n, n, 1));
        boot_medians(i) = median(boot_sample);
    end

    ci = prctile(boot_medians, [100*alpha/2, 100*(1-alpha/2)]);
end

function ci = stratified_bootstrap_ci(data_loc1, data_loc2, n_boot, alpha)
    % Compute bootstrap CI for the average of location-specific medians.
    % This stratified approach resamples within each location independently,
    % computes medians for each location, then averages them. The resulting
    % CI reflects within-location variability, not across-location variability.
    %
    % INPUTS:
    %   data_loc1 - Data vector for location 1
    %   data_loc2 - Data vector for location 2
    %   n_boot    - Number of bootstrap iterations
    %   alpha     - Significance level (e.g., 0.05 for 95% CI)
    %
    % OUTPUT:
    %   ci - [lower, upper] confidence interval

    n1 = length(data_loc1);
    n2 = length(data_loc2);

    boot_avg_medians = nan(n_boot, 1);

    for i = 1:n_boot
        % Resample within each location
        boot_sample1 = data_loc1(randi(n1, n1, 1));
        boot_sample2 = data_loc2(randi(n2, n2, 1));

        % Compute location-specific medians and average them
        med1 = median(boot_sample1);
        med2 = median(boot_sample2);
        boot_avg_medians(i) = (med1 + med2) / 2;
    end

    ci = prctile(boot_avg_medians, [100*alpha/2, 100*(1-alpha/2)]);
end

% --- Condition Filter Functions ---
% Each function takes a trialInfo table and returns a logical column vector

function mask = filter_high_salience(ti)
    % High salience bullseye trials
    mask = (ti.salience == 1);
    mask = mask(:);
end

function mask = filter_low_salience(ti)
    % Low salience bullseye trials
    mask = (ti.salience == 2);
    mask = mask(:);
end

function mask = filter_face(ti)
    % Face image trials
    mask = (ti.stimType == 1);
    mask = mask(:);
end

function mask = filter_nonface(ti)
    % Non-face image trials
    mask = (ti.stimType == 2);
    mask = mask(:);
end

function mask = filter_hirwd_loprob(ti)
    % High reward, low probability trials
    is_high_rwd = (ti.rewardDuration > 200);
    is_low_prob = (ti.blockNumber == 1 & (ti.targetLocIdx == 2 | ti.targetLocIdx == 3 | ti.targetLocIdx == 4)) | ...
                  (ti.blockNumber == 2 & (ti.targetLocIdx == 1 | ti.targetLocIdx == 2 | ti.targetLocIdx == 4));
    mask = is_high_rwd & is_low_prob;
    mask = mask(:);
end

function mask = filter_hirwd_hiprob(ti)
    % High reward, high probability trials
    is_high_rwd = (ti.rewardDuration > 200);
    is_high_prob = (ti.blockNumber == 1 & ti.targetLocIdx == 1) | ...
                   (ti.blockNumber == 2 & ti.targetLocIdx == 3);
    mask = is_high_rwd & is_high_prob;
    mask = mask(:);
end

function mask = filter_lorwd_hiprob(ti)
    % Low reward, high probability trials
    is_low_rwd = (ti.rewardDuration < 200);
    is_high_prob = (ti.blockNumber == 1 & ti.targetLocIdx == 1) | ...
                   (ti.blockNumber == 2 & ti.targetLocIdx == 3);
    mask = is_low_rwd & is_high_prob;
    mask = mask(:);
end

function mask = filter_lorwd_loprob(ti)
    % Low reward, low probability trials
    is_low_rwd = (ti.rewardDuration < 200);
    is_low_prob = (ti.blockNumber == 1 & (ti.targetLocIdx == 2 | ti.targetLocIdx == 3 | ti.targetLocIdx == 4)) | ...
                  (ti.blockNumber == 2 & (ti.targetLocIdx == 1 | ti.targetLocIdx == 2 | ti.targetLocIdx == 4));
    mask = is_low_rwd & is_low_prob;
    mask = mask(:);
end

% --- Main Effect Filter Functions (for Figure 3) ---

function mask = filter_high_reward(ti)
    % All high reward trials (regardless of probability)
    mask = (ti.rewardDuration > 200);
    mask = mask(:);
end

function mask = filter_low_reward(ti)
    % All low reward trials (regardless of probability)
    mask = (ti.rewardDuration < 200);
    mask = mask(:);
end

function mask = filter_high_probability(ti)
    % All high probability trials (regardless of reward)
    % High probability: Block 1 -> Location 1, Block 2 -> Location 3
    mask = (ti.blockNumber == 1 & ti.targetLocIdx == 1) | ...
           (ti.blockNumber == 2 & ti.targetLocIdx == 3);
    mask = mask(:);
end

function mask = filter_low_probability(ti)
    % All low probability trials (regardless of reward)
    % Low probability: Block 1 -> Locations 2,3,4; Block 2 -> Locations 1,2,4
    mask = (ti.blockNumber == 1 & (ti.targetLocIdx == 2 | ti.targetLocIdx == 3 | ti.targetLocIdx == 4)) | ...
           (ti.blockNumber == 2 & (ti.targetLocIdx == 1 | ti.targetLocIdx == 2 | ti.targetLocIdx == 4));
    mask = mask(:);
end
