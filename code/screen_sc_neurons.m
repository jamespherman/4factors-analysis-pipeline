function [selected_neurons, sig_epoch_comparison, ...
    scSide, neuron_class, session_data, ...
    neuron_class_v2] = ...
    screen_sc_neurons(session_data, project_root, ...
    condition_defs, session_name)
% screen_sc_neurons - Classifies SC neurons into Classic, Interneuron, or Excluded.
%
% This function implements a three-tier classification system for SC neurons:
%
% Significance and rates use successful gSac_4factors trials inside
% each neuron's stable period. Modulation = paired t-test of epoch vs
% baseline rate per location x epoch, Bonferroni-corrected (see
% define_task_conditions, SC inclusion gate).
%
% CLASSIC SC NEURON (all must be true):
%   (a) Significant increase in visual OR delay OR saccade epoch vs baseline
%       for at least one contralateral location
%   (b) Peak evoked rate >= min_evoked_rate_sc (5 sp/s): trial-averaged
%       visual, delay or saccade rate at the best location
%   (c) Every significant ipsilateral activation (any epoch) is
%       reliably less than half of contra: the best contralateral
%       location's per-trial response in that epoch is significantly
%       larger than twice the ipsilateral one (Welch, one-sided, 0.05)
%   (d) NO significant suppression (any location) in the visual or
%       saccade epoch; delay suppression (mostly ramp-down) is exempt
%   (Rule decided 2026-10-28, shared SC spec section 6; replaces the
%   opposing-location and contralateral-suppression criteria.)
%
% PUTATIVE INTERNEURON (non-classic):
%   - Peak evoked rate >= min_evoked_rate_sc
%   - Task-modulated (activation or suppression somewhere)
%   - Fails (a), (c) or (d)
%
% EXCLUDED:
%   - Peak evoked rate < min_evoked_rate_sc
%   - Not task-modulated
%
% Mean FR and sparsity are still computed (stable period) but no
% longer gate inclusion (changed 2026-10-16).
%
% INPUTS:
%   session_data   - A struct containing session-specific data, conforming to the
%                    `session_data_dictionary.md`.
%   project_root   - The root path of the project directory.
%   condition_defs - (Optional) The analysis plan struct from define_task_conditions.
%                    If condition_defs.neuron_inclusion.use_strict_screening is
%                    false, all neurons are selected (bypassing screening criteria).
%   session_name   - (Optional) Session identifier for figure filenames.
%
% OUTPUT:
%   selected_neurons     - A logical vector (nClusters x 1) where true indicates
%                          a neuron classified as 'classic' or 'interneuron'.
%   sig_epoch_comparison - A logical matrix (nClusters x 3) indicating
%                          significant firing rate changes between epochs.
%   scSide               - A string ('right' or 'left') indicating the determined
%                          recorded SC side.
%   neuron_class         - A cell array (nClusters x 1) with classification:
%                          'classic', 'interneuron', or 'excluded'.
%   session_data         - Updated session_data with cached_friedman_results
%                          for faster subsequent runs.
%   neuron_class_v2      - (Optional) Cell array (nClusters x 1) with
%                          data-driven cluster labels from v2
%                          classification. Empty if clustering not run.
%

% Check if screening should be bypassed
use_strict_screening = true; % Default to strict screening for backward compatibility
if nargin >= 3 && ~isempty(condition_defs) && ...
        isfield(condition_defs, 'neuron_inclusion') && ...
        isfield(condition_defs.neuron_inclusion, 'use_strict_screening')
    use_strict_screening = condition_defs.neuron_inclusion.use_strict_screening;
end

% Determine session_id for figure filenames
if nargin >= 4 && ~isempty(session_name)
    session_id = session_name;
elseif isfield(session_data, 'metadata') && isfield(session_data.metadata, 'unique_id')
    session_id = session_data.metadata.unique_id;
else
    session_id = 'unknown_session';
end

fprintf('screen_sc_neurons: Identifying task-modulated neurons...\n');

% Define output directory and filename
output_dir = fullfile(project_root, 'figures');

% --- Setup and Data Extraction ---
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info.cluster_id);

all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

trialInfo = session_data.trialInfo;
eventTimes = session_data.eventTimes;
codes = initCodes();

% --- Load Stability Info (for stable period calculations) ---
% stability_info contains per-neuron stable recording periods identified by
% detect_stable_trial_period.m. If available, we use it to calculate metrics
% only over each neuron's stable period.
use_stable_periods = false;
stability_info = [];
if isfield(session_data, 'stability_info') && ~isempty(session_data.stability_info)
    stability_info = session_data.stability_info;
    use_stable_periods = true;
    fprintf('screen_sc_neurons: Using stability_info for per-neuron stable periods.\n');
else
    fprintf('screen_sc_neurons: No stability_info found, using full session for metrics.\n');
end

% Initialize outputs for early return
selected_neurons = false(nClusters, 1);
sig_epoch_comparison = false(nClusters, 3);
scSide = 'unknown';
neuron_class = repmat({'excluded'}, nClusters, 1);
neuron_class_v2 = repmat({'excluded'}, ...
    nClusters, 1);

% Classification thresholds
% Mean FR and sparsity: reported only, no longer gate inclusion
MIN_FR_THRESHOLD = 3.0;  % sp/s - minimum mean firing rate (changed from 3.5)
MAX_SPARSE_FRACTION = 0.80;  % Maximum fraction of empty 100ms bins (changed from 0.70)

% Inclusion gate (define_task_conditions.neuron_inclusion)
MIN_EVOKED_RATE = 5.0;
MOD_ALPHA = 0.05;
if nargin >= 3 && ~isempty(condition_defs) && ...
        isfield(condition_defs, 'neuron_inclusion')
    ni = condition_defs.neuron_inclusion;
    if isfield(ni, 'min_evoked_rate_sc')
        MIN_EVOKED_RATE = ni.min_evoked_rate_sc;
    end
    if isfield(ni, 'sc_modulation_alpha')
        MOD_ALPHA = ni.sc_modulation_alpha;
    end
end
% Classic vs interneuron: one-sided Welch alpha of the test that the
% contralateral response exceeds twice an ipsilateral activation
IPSI_ALPHA = 0.05;
% Cached test results are reused only if computed by this method
TEST_METHOD = 'paired_ttest_stable4f_v1';

if nClusters == 0
    fprintf('WARNING in screen_sc_neurons: No clusters found.\n');
    return;
end

%% 2. Trial Identification
% Find all valid memory-guided saccade trials from both gSac_jph and
% gSac_4factors tasks. These are trials with a recorded target
% re-illumination and a reward, indicating successful completion.

% Identify gSac_jph memory-guided saccade trials
gSac_jph_memSac_trials = find(trialInfo.taskCode == ...
    codes.uniqueTaskCode_gSac_jph & ...
    eventTimes.targetReillum > 0 & ...
    eventTimes.pdsReward > 0);
fprintf('Found %d valid gSac_jph memory-guided trials.\n', ...
    length(gSac_jph_memSac_trials));

% Identify gSac_4factors memory-guided saccade trials
gSac_4factors_memSac_trials = find(trialInfo.taskCode == ...
    codes.uniqueTaskCode_gSac_4factors & ...
    ~isnan(eventTimes.targetReillum) & ...
    eventTimes.pdsReward > 0);
fprintf('Found %d valid gSac_4factors memory-guided trials.\n', ...
    length(gSac_4factors_memSac_trials));

% Combine all trials for a unified firing rate calculation
all_memSac_trials = union(gSac_jph_memSac_trials, ...
    gSac_4factors_memSac_trials);

if isempty(all_memSac_trials)
    fprintf(['WARNING in screen_sc_neurons: No suitable memory-guided ' ...
        'trials found in either gSac_jph or gSac_4factors tasks. ' ...
        'Skipping.\n']);
    return;
end

nMemSacTrials = length(all_memSac_trials);

% If use_strict_screening is false, include all neurons and skip detailed analysis
if ~use_strict_screening
    fprintf('screen_sc_neurons: use_strict_screening=false, including all %d neurons.\n', nClusters);
    selected_neurons = true(nClusters, 1);
    sig_epoch_comparison = false(nClusters, 3); % No significance testing performed
    neuron_class = repmat({'classic'}, nClusters, 1); % Label all as classic when bypassing

    % Determine scSide using the primary method (gSac_jph trials) only
    if length(gSac_jph_memSac_trials) > 5
        thetas_jph = trialInfo.targetTheta(gSac_jph_memSac_trials) / 10;
        left_vf_trials = sum(thetas_jph > 90 & thetas_jph < 270);
        right_vf_trials = sum(thetas_jph < 90 | thetas_jph > 270);
        if left_vf_trials > right_vf_trials
            scSide = 'right';
        else
            scSide = 'left';
        end
        fprintf('Determined SC Side: %s (from gSac_jph trials).\n', scSide);
    else
        % Fallback: use grid_hole from metadata if available
        if isfield(session_data, 'metadata') && isfield(session_data.metadata, 'grid_hole')
            grid_hole_str = session_data.metadata.grid_hole;
            coords = sscanf(grid_hole_str, '(%f, %f)');
            if ~isempty(coords)
                if coords(1) < 0
                    scSide = 'left';
                elseif coords(1) > 0
                    scSide = 'right';
                else
                    scSide = 'unknown';
                end
                fprintf('Determined SC Side: %s (from grid_hole metadata).\n', scSide);
            else
                scSide = 'unknown';
                fprintf('Could not determine SC Side.\n');
            end
        else
            scSide = 'unknown';
            fprintf('Could not determine SC Side (no gSac_jph trials or metadata).\n');
        end
    end

    fprintf('Finished screening (all neurons included). Found %d neurons.\n', nClusters);
    return;
end

%% 3. Vectorized Firing Rate Calculation
% Epoch definitions: {event_name, start_offset, end_offset, duration}
epochs = {
    'targetOn',     -0.075, 0.025,  0.1;   % 1. Baseline
    'targetOn',     0.05,   0.2,    0.15;  % 2. Visual
    'fixOff',       -0.15,  0.05,   0.2;   % 3. Delay
    'saccadeOnset', -0.025, 0.05,   0.075  % 4. Saccade
    };
nEpochs = size(epochs, 1);
epoch_frs = nan(nClusters, nEpochs, nMemSacTrials);

% This new neuron-centric approach iterates through each neuron and calculates
% its firing rate for all trials and all epochs in a vectorized manner,
% which is more efficient than iterating through each trial.
for i_cluster = 1:nClusters
    % Get all spike times for the current cluster
    spike_times = all_spike_times(all_spike_clusters == ...
        cluster_ids(i_cluster));

    if isempty(spike_times)
        % If no spikes, FR is 0 for all epochs and trials for this cluster
        epoch_frs(i_cluster, :, :) = 0;
        continue;
    end

    for i_epoch = 1:nEpochs
        event_name = epochs{i_epoch, 1};
        win_dur    = epochs{i_epoch, 4};

        % Get event times for all relevant trials for the current epoch
        epoch_event_times = eventTimes.(event_name)(all_memSac_trials);

        % Create a [nMemSacTrials x 2] matrix of time windows
        time_windows = [epoch_event_times + epochs{i_epoch, 2}, ...
            epoch_event_times + epochs{i_epoch, 3}];

        % Find trials where the event time is NaN and keep track of them
        valid_trials_mask = ~isnan(epoch_event_times);

        % Filter out trials with NaN event times
        valid_time_windows = time_windows(valid_trials_mask, :);

        if isempty(valid_time_windows)
            continue; % No valid trials for this epoch, NaNs will remain
        end

        % Reshape the time_windows matrix into a single row vector of edges
        % for histcounts: [start1, end1, start2, end2, ...]
        edges = reshape(valid_time_windows', 1, []);

        try
            % Get the counts for all bins (both inside and outside the windows)
            all_counts = histcounts(spike_times, edges);
        catch me
            keyboard
        end

        % The counts within our desired windows are the odd-indexed elements
        % (1st, 3rd, 5th, etc.) of the histcounts output.
        spike_counts_in_windows = all_counts(1:2:end);

        try
            % Create a temporary array to store results for the current epoch
            temp_frs = nan(1, nMemSacTrials);
            temp_frs(valid_trials_mask) = spike_counts_in_windows / win_dur;
        catch me
            keyboard
        end

        % Place the calculated firing rates into the master matrix
        epoch_frs(i_cluster, i_epoch, :) = temp_frs;
    end
end

%% 4. Calculate Session-Wide Metrics (Mean FR and Sparseness)
% These metrics are needed for the three-tier classification system.
% If stability_info is available, calculate metrics only over each neuron's
% stable recording period. Otherwise, use the full session.

% Determine full session time range (for fallback)
session_start = min(all_spike_times);
session_end = max(all_spike_times);
session_duration = session_end - session_start;

% Calculate mean FR and sparseness for each neuron
mean_fr_session = nan(nClusters, 1);
sparse_fraction = nan(nClusters, 1);
bin_size = 0.1; % 100ms bins for sparseness calculation

n_using_stable = 0;
for i_cluster = 1:nClusters
    % Get all spike times for this cluster
    spike_times = all_spike_times(all_spike_clusters == cluster_ids(i_cluster));

    % Determine time range for this neuron
    % Use stable period if available, otherwise full session
    if use_stable_periods && i_cluster <= length(stability_info) && ...
            isfield(stability_info(i_cluster), 'is_stable') && ...
            stability_info(i_cluster).is_stable
        % Get stable trial range and convert to times
        stable_trial_start = stability_info(i_cluster).stable_trial_start;
        stable_trial_end = stability_info(i_cluster).stable_trial_end;

        if ~isnan(stable_trial_start) && ~isnan(stable_trial_end) && ...
                stable_trial_start >= 1 && stable_trial_end >= stable_trial_start && ...
                stable_trial_end <= length(eventTimes.trialBegin)
            % Convert trial indices to times
            neuron_start = eventTimes.trialBegin(stable_trial_start);
            neuron_end = eventTimes.trialEnd(stable_trial_end);
            n_using_stable = n_using_stable + 1;
        else
            % Invalid stable period, use full session
            neuron_start = session_start;
            neuron_end = session_end;
        end
    else
        % No stable period info, use full session
        neuron_start = session_start;
        neuron_end = session_end;
    end

    neuron_duration = neuron_end - neuron_start;

    % Filter spike times to the neuron's active period
    spike_times_filtered = spike_times(spike_times >= neuron_start & spike_times <= neuron_end);

    % Mean firing rate across the neuron's active period
    if neuron_duration > 0
        mean_fr_session(i_cluster) = length(spike_times_filtered) / neuron_duration;
    else
        mean_fr_session(i_cluster) = 0;
    end

    % Sparseness: fraction of 100ms bins that are empty within active period
    if ~isempty(spike_times_filtered) && neuron_duration > 0
        n_bins = ceil(neuron_duration / bin_size);
        bin_edges = neuron_start + (0:n_bins) * bin_size;
        spike_counts = histcounts(spike_times_filtered, bin_edges);
        empty_bins = sum(spike_counts == 0);
        sparse_fraction(i_cluster) = empty_bins / n_bins;
    else
        sparse_fraction(i_cluster) = 1; % No spikes = 100% sparse
    end
end

if use_stable_periods
    fprintf('Metrics calculated using stable periods for %d/%d neurons.\n', n_using_stable, nClusters);
end
fprintf('Mean FR range: %.2f - %.2f sp/s\n', min(mean_fr_session), max(mean_fr_session));
fprintf('Sparseness range: %.1f%% - %.1f%% empty bins\n', ...
    min(sparse_fraction)*100, max(sparse_fraction)*100);

%% 5. Hierarchical scSide Determination (must run before classification)
% Determine the recorded side of the SC. The primary method uses gSac_jph
% trials, as the experimenter-placed targets are considered ground truth.
% If insufficient gSac_jph trials exist, a fallback method uses
% gSac_4factors data to compare population-level visual responses.

if length(gSac_jph_memSac_trials) > 5
    % Primary method: Use gSac_jph trials
    % Because these trials are placed in the contralateral field by the
    % experimenter, we can determine scSide based on target locations.
    thetas_jph = trialInfo.targetTheta(gSac_jph_memSac_trials) / 10;
    left_vf_trials = sum(thetas_jph > 90 & thetas_jph < 270);
    right_vf_trials = sum(thetas_jph < 90 | thetas_jph > 270);

    if left_vf_trials > right_vf_trials
        scSide = 'right'; % Right SC records from the left visual field
    else
        scSide = 'left';  % Left SC records from the right visual field
    end
    fprintf(['Determined SC Side: %s based on contralateral target ' ...
        'placement in gSac_jph task.\n'], scSide);
else
    % Fallback method: Use gSac_4factors trials
    % Compare population average visual response for left vs. right targets.

    % Create a logical mask for which of the combined trials belong to the
    % gSac_4factors task.
    is_4factors_trial = ismember(all_memSac_trials, gSac_4factors_memSac_trials);

    % Get target angles for all trials in the combined list.
    thetas_all = trialInfo.targetTheta(all_memSac_trials) / 10;

    % Create masks for left and right hemifield trials, but only apply them
    % to the gSac_4factors trials.
    left_trials_mask = is_4factors_trial & (thetas_all > 90 & thetas_all < 270);
    right_trials_mask = is_4factors_trial & (thetas_all < 90 | thetas_all > 270);

    % Calculate avg visual FR for left vs. right trials across all neurons.
    mean_vis_fr_left = mean(epoch_frs(:, 2, left_trials_mask), 3, 'omitnan');
    mean_vis_fr_right = mean(epoch_frs(:, 2, right_trials_mask), 3, 'omitnan');

    % Compare the mean across the entire population to determine side.
    if mean(mean_vis_fr_left, 'omitnan') > mean(mean_vis_fr_right, 'omitnan')
        scSide = 'right'; % Right SC prefers left visual field
    else
        scSide = 'left'; % Left SC prefers right visual field
    end
    fprintf(['Insufficient gSac_jph trials. Determined SC Side: %s by ' ...
        'comparing population visual responses in gSac_4factors.\n'], scSide);
end

%% 6. Three-Tier Neuron Classification
% Classify each neuron as 'classic', 'interneuron', or 'excluded' based on:
%
% CLASSIC SC NEURON (all must be true):
%   (a) Significant increase vs baseline for >=1 contralateral location
%   (b) Mean FR > 3.0 sp/s
%   (c) NOT significantly activated at opposing locations (180 deg apart)
%   (d) Not sparse (<=80% of 100ms bins empty)
%
% PUTATIVE INTERNEURON:
%   - Mean FR > 3.0 sp/s
%   - Task-modulated (some significant effect somewhere)
%   - Fails to meet "classic" criteria
%
% EXCLUDED:
%   - Mean FR <= 3.0 sp/s OR not task-modulated OR sparse

% Define contralateral locations based on scSide
% Right SC responds to left visual field: 90 < theta < 270
% Left SC responds to right visual field: theta < 90 OR theta > 270
if strcmp(scSide, 'right')
    is_contralateral = @(theta_deg) theta_deg > 90 & theta_deg < 270;
elseif strcmp(scSide, 'left')
    is_contralateral = @(theta_deg) theta_deg < 90 | theta_deg > 270;
else
    % Unknown side - treat all locations as contralateral
    is_contralateral = @(theta_deg) true(size(theta_deg));
    fprintf('WARNING: scSide unknown, treating all locations as contralateral.\n');
end

% Find logical indices for each task within the combined trial array
is_jph_trial = ismember(all_memSac_trials, gSac_jph_memSac_trials);
is_4factors_trial = ismember(all_memSac_trials, gSac_4factors_memSac_trials);

% Get the unique target locations for the 4factors task (in degrees)
unique_locations_4factors = unique(trialInfo.targetTheta(gSac_4factors_memSac_trials));
unique_locations_deg = unique_locations_4factors / 10;
n_locations = length(unique_locations_4factors);

% Get target thetas for all combined trials
thetas_all = trialInfo.targetTheta(all_memSac_trials);

% Build trial groups - one per unique location
trial_groups = cell(1, n_locations);
location_angles = nan(1, n_locations);
for i_loc = 1:n_locations
    loc = unique_locations_4factors(i_loc);
    trial_groups{i_loc} = (thetas_all == loc) & is_4factors_trial;
    location_angles(i_loc) = loc / 10; % degrees
end

% Also add gSac_jph trials as a separate group if available
has_jph = any(is_jph_trial);
if has_jph
    trial_groups{end+1} = is_jph_trial;
    % For jph trials, get the most common target angle
    jph_thetas = trialInfo.targetTheta(gSac_jph_memSac_trials) / 10;
    location_angles(end+1) = mode(jph_thetas);
end

n_groups = length(trial_groups);

% Initialize storage for per-location significance results
% sig_by_location: nClusters x n_groups x 3 (vis, delay, sac)
sig_by_location = false(nClusters, n_groups, 3);
% suppression_by_location: tracks significant DECREASES vs baseline
suppression_by_location = false(nClusters, n_groups, 3);
group_mean_frs = cell(1, n_groups);
for i_group = 1:n_groups
    group_mean_frs{i_group} = nan(nClusters, nEpochs);
end
% Peak evoked rate: best 4factors location, visual/delay/saccade epoch
peak_evoked_rate = zeros(nClusters, 1);
n_4f_groups = n_locations;

% Per-neuron trial mask over all_memSac_trials: 4factors trials inside
% the neuron's stable period (all 4factors trials if it has none).
% Stable periods are scoped to gSac_4factors, so gSac_jph trials
% outside them drop out of the tests.
stable_ranges = nan(nClusters, 2);
neuron_trial_ok = repmat(is_4factors_trial(:)', nClusters, 1);
if use_stable_periods
    for i_cluster = 1:min(nClusters, length(stability_info))
        si_c = stability_info(i_cluster);
        if isfield(si_c, 'is_stable') && si_c.is_stable && ...
                ~isnan(si_c.stable_trial_start)
            stable_ranges(i_cluster, :) = ...
                [si_c.stable_trial_start, si_c.stable_trial_end];
            neuron_trial_ok(i_cluster, :) = ...
                all_memSac_trials(:)' >= si_c.stable_trial_start & ...
                all_memSac_trials(:)' <= si_c.stable_trial_end;
        end
    end
end

%% Check for cached test results
% Reused only if computed by the same method on the same stable periods
use_cached_friedman = false;
if isfield(session_data, 'cached_friedman_results') && ...
        ~isempty(session_data.cached_friedman_results)
    cached = session_data.cached_friedman_results;
    if isequal(cached.cluster_ids(:), cluster_ids(:)) && ...
            size(cached.sig_by_location, 2) == n_groups && ...
            isfield(cached, 'suppression_by_location') && ...
            isfield(cached, 'method') && ...
            strcmp(cached.method, TEST_METHOD) && ...
            isequaln(cached.stable_ranges, stable_ranges)
        use_cached_friedman = true;
        fprintf('screen_sc_neurons: Using cached test results.\n');
        sig_by_location = cached.sig_by_location;
        suppression_by_location = cached.suppression_by_location;
        group_mean_frs = cached.group_mean_frs;
        peak_evoked_rate = cached.peak_evoked_rate;
    else
        fprintf(['screen_sc_neurons: Cached test results stale ' ...
            '(method, clusters or stable periods changed), ' ...
            'recomputing...\n']);
    end
end

if ~use_cached_friedman
% Bonferroni over every location x epoch test
alpha_corr = MOD_ALPHA / (3 * n_groups);
for i_cluster = 1:nClusters
    for i_group = 1:n_groups
        trial_mask = trial_groups{i_group}(:)' & ...
            neuron_trial_ok(i_cluster, :);

        if sum(trial_mask) < 2
            continue;
        end

        % Extract firing rates for the current neuron and trial group
        neuron_frs_all_trials = squeeze(epoch_frs(i_cluster, :, trial_mask))';
        neuron_frs = neuron_frs_all_trials(~any(isnan(neuron_frs_all_trials), 2), :);

        if size(neuron_frs, 1) < 5
            continue;
        end

        mean_fr_this_group = mean(neuron_frs, 1, 'omitnan');
        group_mean_frs{i_group}(i_cluster, :) = mean_fr_this_group;
        if i_group <= n_4f_groups
            peak_evoked_rate(i_cluster) = max(peak_evoked_rate( ...
                i_cluster), max(mean_fr_this_group(2:4)));
        end

        % Paired t-test of each epoch vs baseline (Vis, Del, Sac). The
        % epochs differ in length, so a test of the mean is used; rank
        % tests flag unmodulated low-rate neurons in that case.
        for i_comp = 1:3
            d = neuron_frs(:, i_comp + 1) - neuron_frs(:, 1);
            if all(d == 0), continue; end
            [~, p] = ttest(d);
            if p < alpha_corr
                if mean(d) > 0
                    % Significant INCREASE vs baseline
                    sig_by_location(i_cluster, i_group, i_comp) = true;
                else
                    % Significant DECREASE (suppression) vs baseline
                    suppression_by_location(i_cluster, i_group, ...
                        i_comp) = true;
                end
            end
        end
    end
end

    % Cache results in session_data for future runs
    session_data.cached_friedman_results = struct(...
        'cluster_ids', cluster_ids, ...
        'sig_by_location', sig_by_location, ...
        'suppression_by_location', suppression_by_location, ...
        'group_mean_frs', {group_mean_frs}, ...
        'location_angles', location_angles, ...
        'epoch_frs', epoch_frs, ...
        'peak_evoked_rate', peak_evoked_rate, ...
        'method', TEST_METHOD, ...
        'stable_ranges', stable_ranges, ...
        'computation_date', datestr(now));
    fprintf('screen_sc_neurons: Test results cached for future runs.\n');
end  % end of if ~use_cached_friedman

% Now classify each neuron
for i_cluster = 1:nClusters
    % Evoked-rate criterion (mean FR and sparsity no longer gate)
    has_evoked = peak_evoked_rate(i_cluster) >= MIN_EVOKED_RATE;

    % Determine if neuron is task-modulated at any location (activation or suppression)
    is_task_modulated = any(sig_by_location(i_cluster, :, :), 'all') || ...
                        any(suppression_by_location(i_cluster, :, :), 'all');

    % Check for significant activation at contralateral locations
    % Also track activated locations for opposing activation check
    has_contralateral_activation = false;
    activated_locations = []; % track which locations show ANY activation
    activated_locations_vis_sac = []; % Only Visual (1) and Saccade (3) epochs

    for i_group = 1:n_groups
        theta = location_angles(i_group);

        % Check if this location has any activation (any epoch)
        if any(sig_by_location(i_cluster, i_group, :))
            activated_locations(end+1) = theta;
            if is_contralateral(theta)
                has_contralateral_activation = true;
            end
        end

        % Track Visual (epoch 1) or Saccade (epoch 3) activation separately
        % for opposing activation check (NOT delay epoch 2)
        if sig_by_location(i_cluster, i_group, 1) || sig_by_location(i_cluster, i_group, 3)
            activated_locations_vis_sac(end+1) = theta;
        end
    end

    % Ipsilateral activation, rule C with a half-of-contra null (shared
    % SC spec section 6; decided 2026-10-28). An ipsilateral activation
    % (any epoch) makes the neuron non-classic unless, in the same
    % epoch, the best contralateral location's per-trial response
    % (epoch - baseline) is significantly larger than TWICE the
    % ipsilateral one (Welch, one-sided, IPSI_ALPHA), i.e. the ipsi
    % response is reliably less than half of contra. Plain significance
    % (rule A) depends on trial count: with ~300 trials per location,
    % responses at a fraction of contra become significant.
    has_ipsi_activation = false;
    resp = squeeze(epoch_frs(i_cluster, 2:4, :))' - ...
        squeeze(epoch_frs(i_cluster, 1, :));
    con_groups = find(arrayfun(is_contralateral, ...
        location_angles(1:n_4f_groups)));
    for i_group = 1:n_4f_groups
        if is_contralateral(location_angles(i_group))
            continue;
        end
        for i_ep = find(squeeze(sig_by_location(i_cluster, i_group, ...
                :)))'
            ip = resp(trial_groups{i_group}(:) & ...
                neuron_trial_ok(i_cluster, :)', i_ep);
            ip = ip(~isnan(ip));
            best_p = 1;
            for i_con = con_groups
                cp = resp(trial_groups{i_con}(:) & ...
                    neuron_trial_ok(i_cluster, :)', i_ep);
                cp = cp(~isnan(cp));
                if numel(cp) < 3 || numel(ip) < 3, continue; end
                [~, p] = ttest2(cp, 2 * ip, 'Vartype', 'unequal', ...
                    'Tail', 'right');
                best_p = min(best_p, p);
            end
            if best_p >= IPSI_ALPHA
                has_ipsi_activation = true;
            end
        end
    end

    % Suppression (any location) in the visual or saccade epoch. Delay
    % suppression is exempt: against the pre-target baseline it is
    % mostly a delay ramp-down, which James counts as classic (Newton
    % precedent, 2026-10-21).
    has_suppression = any(suppression_by_location(i_cluster, ...
        1:n_4f_groups, [1 3]), 'all');

    % Apply classification logic
    is_interneuron_pattern = has_ipsi_activation || has_suppression;

    if ~has_evoked || ~is_task_modulated
        % EXCLUDED: weak evoked response or not task-modulated
        neuron_class{i_cluster} = 'excluded';
        selected_neurons(i_cluster) = false;
    elseif has_contralateral_activation && ~is_interneuron_pattern
        % CLASSIC: contralateral activation, no counted ipsilateral
        % activation, no visual/saccade suppression
        neuron_class{i_cluster} = 'classic';
        selected_neurons(i_cluster) = true;
    else
        % INTERNEURON (non-classic): evoked response, task-modulated,
        % but ipsilateral activation, suppression, or no contralateral
        % activation
        neuron_class{i_cluster} = 'interneuron';
        selected_neurons(i_cluster) = true;
    end

    % Store canonical significance profile (from first significant location)
    for i_group = 1:n_groups
        if any(sig_by_location(i_cluster, i_group, :))
            sig_epoch_comparison(i_cluster, :) = squeeze(sig_by_location(i_cluster, i_group, :))';
            break;
        end
    end
end

%% Apply Manual Classification Overrides
% Some neurons are misclassified by the automated criteria due to edge cases
% (e.g., weak but statistically significant ipsilateral responses triggering
% opposing activation). These overrides are applied after visual inspection
% of individual neuron summary PDFs.
%
% Format: {session_name, cluster_id, forced_class, reason}
% (N14 cl 6 was forced to classic here; James judged it non-classic on
% 2026-10-23, which the 2026-10-28 rule gives without an override.)
manual_overrides = cell(0, 4);

for i_override = 1:size(manual_overrides, 1)
    override_session = manual_overrides{i_override, 1};
    override_cluster = manual_overrides{i_override, 2};
    override_class = manual_overrides{i_override, 3};
    override_reason = manual_overrides{i_override, 4};

    % Check if this override applies to current session
    if strcmp(session_id, override_session)
        % Find the cluster index
        override_idx = find(cluster_ids == override_cluster);
        if ~isempty(override_idx)
            old_class = neuron_class{override_idx};
            neuron_class{override_idx} = override_class;
            % Update selected_neurons based on new class
            if strcmp(override_class, 'classic') || strcmp(override_class, 'interneuron')
                selected_neurons(override_idx) = true;
            else
                selected_neurons(override_idx) = false;
            end
            fprintf('MANUAL OVERRIDE: %s cluster %d: %s -> %s (%s)\n', ...
                override_session, override_cluster, old_class, override_class, override_reason);
        end
    end
end

% Count classifications
n_classic = sum(strcmp(neuron_class, 'classic'));
n_interneuron = sum(strcmp(neuron_class, 'interneuron'));
n_excluded = sum(strcmp(neuron_class, 'excluded'));

fprintf('\n=== SC Neuron Classification Results ===\n');
fprintf('Classic SC neurons:    %d\n', n_classic);
fprintf('Putative interneurons: %d\n', n_interneuron);
fprintf('Excluded:              %d\n', n_excluded);
fprintf('Total selected:        %d / %d\n', nnz(selected_neurons), nClusters);

% --- Generate summary figure ---
if n_groups > 0
    % Calculate global max firing rate for consistent color scaling
    global_max_fr = 0;
    for i_group = 1:n_groups
        max_in_group = max(group_mean_frs{i_group}, [], 'all', 'omitnan');
        if max_in_group > global_max_fr
            global_max_fr = max_in_group;
        end
    end
    if global_max_fr == 0; global_max_fr = 1; end

    fig = figure('Color', 'w', 'Position', [100, 100, 350 * n_groups, 700]);

    fr_x_labels = {'Base', 'Vis', 'Delay', 'Sac'};
    sig_x_labels = {'Vis', 'Delay', 'Sac'};

    for i_group = 1:n_groups
        % Subplot for Firing Rates
        plot_idx_fr = (i_group - 1) * 2 + 1;
        mySubPlot([1, n_groups * 2, plot_idx_fr], ...
            'Width', 0.92, 'LeftMargin', 0.05);
        imagesc(group_mean_frs{i_group});
        clim([0, global_max_fr]);
        colormap(gca, flipud(hot));
        set(gca, 'XTick', 1:4, 'XTickLabel', fr_x_labels);
        if plot_idx_fr == 1
            ylabel('Cluster ID');
        else
            set(gca, 'YTickLabel', []);
        end

        % Generate title
        theta = location_angles(i_group);
        contra_str = '';
        if is_contralateral(theta)
            contra_str = ' (contra)';
        end
        if has_jph && i_group == n_groups
            title_str = sprintf('gSac_jph: %d%s%s', theta, char(176), contra_str);
        else
            title_str = sprintf('%d%s%s', theta, char(176), contra_str);
        end
        title(title_str);

        % Subplot for Significance
        plot_idx_sig = (i_group - 1) * 2 + 2;
        mySubPlot([1, n_groups * 2, plot_idx_sig], 'Width', 0.92, ...
            'LeftMargin', 0.05);
        imagesc(squeeze(sig_by_location(:, i_group, :)));
        colormap(gca, flipud(bone));
        set(gca, 'XTick', 1:3, 'XTickLabel', sig_x_labels, 'YTickLabel', []);
    end

    % Save the figure
    figFileName = fullfile(output_dir, [session_id, '_sc_epoch_frs.pdf']);
    pdfSave(figFileName, fig.Position(3:4)/72, fig);
    close(fig);
end

fprintf(['Finished screening. Selected %d neurons ' ...
    '(Classic: %d, Interneuron: %d).\n'], ...
    nnz(selected_neurons), n_classic, n_interneuron);

%% Data-Driven Classification (v2)
% Extract features and cluster neurons that pass
% quality gates. This runs alongside v1 and does NOT
% affect v1 outputs (selected_neurons, neuron_class).
if nargout >= 6
    try
        features = extract_sc_response_features( ...
            session_data, scSide);
        n_quality = sum(features.passes_quality);

        if n_quality >= 30
            clustering = cluster_sc_neurons( ...
                features.feature_matrix( ...
                    features.passes_quality, :), ...
                features.feature_names, ...
                features.cluster_ids( ...
                    features.passes_quality));

            % Map cluster assignments back
            quality_idx = find( ...
                features.passes_quality);
            for i_q = 1:length(quality_idx)
                neuron_class_v2{quality_idx(i_q)} ...
                    = clustering.labels{ ...
                    clustering.assignments(i_q)};
            end

            fprintf(['v2 classification: ' ...
                '%d neurons in %d clusters.\n'], ...
                n_quality, clustering.k);
        else
            % Too few neurons for clustering
            fprintf(['v2 classification: only ' ...
                '%d quality neurons (need 30). ' ...
                'Falling back to v1.\n'], ...
                n_quality);
            neuron_class_v2 = neuron_class;
        end
    catch ME
        warning(['v2 classification failed: ' ...
            '%s. Using v1 fallback.'], ...
            ME.message);
        neuron_class_v2 = neuron_class;
    end
end

end
