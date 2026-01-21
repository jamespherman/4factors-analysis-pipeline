function [selected_neurons, sig_epoch_comparison, scSide, neuron_class] = screen_sc_neurons(session_data, project_root, condition_defs)
% screen_sc_neurons - Classifies SC neurons into Classic, Interneuron, or Excluded.
%
% This function implements a three-tier classification system for SC neurons:
%
% CLASSIC SC NEURON (all must be true):
%   (a) Significant increase in visual OR delay OR saccade epoch vs baseline
%       for at least one contralateral location
%   (b) Mean FR > 3.5 sp/s
%   (c) NOT significantly activated at opposing locations (e.g., not both
%       135 deg and 315 deg)
%   (d) Not sparse (<=70% of 100ms bins empty across session)
%
% PUTATIVE INTERNEURON:
%   - Mean FR > 3.5 sp/s
%   - Task-modulated (some significant effect somewhere)
%   - Fails to meet "classic" criteria (bilateral responses, suppression, etc.)
%
% EXCLUDED:
%   - Mean FR <= 3.5 sp/s
%   - Not task-modulated
%   - Sparse (>70% of 100ms bins empty)
%
% INPUTS:
%   session_data   - A struct containing session-specific data, conforming to the
%                    `session_data_dictionary.md`.
%   project_root   - The root path of the project directory.
%   condition_defs - (Optional) The analysis plan struct from define_task_conditions.
%                    If condition_defs.neuron_inclusion.use_strict_screening is
%                    false, all neurons are selected (bypassing screening criteria).
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
%

% Check if screening should be bypassed
use_strict_screening = true; % Default to strict screening for backward compatibility
if nargin >= 3 && ~isempty(condition_defs) && ...
        isfield(condition_defs, 'neuron_inclusion') && ...
        isfield(condition_defs.neuron_inclusion, 'use_strict_screening')
    use_strict_screening = condition_defs.neuron_inclusion.use_strict_screening;
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

% Initialize outputs for early return
selected_neurons = false(nClusters, 1);
sig_epoch_comparison = false(nClusters, 3);
scSide = 'unknown';
neuron_class = repmat({'excluded'}, nClusters, 1);

% Classification thresholds
MIN_FR_THRESHOLD = 3.5;  % sp/s - minimum mean firing rate
MAX_SPARSE_FRACTION = 0.70;  % Maximum fraction of empty 100ms bins

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

% Determine session time range from all spike times
session_start = min(all_spike_times);
session_end = max(all_spike_times);
session_duration = session_end - session_start;

% Calculate mean FR and sparseness for each neuron
mean_fr_session = nan(nClusters, 1);
sparse_fraction = nan(nClusters, 1);

% Create 100ms bins spanning the session
bin_size = 0.1; % 100ms
n_bins = ceil(session_duration / bin_size);
bin_edges = session_start + (0:n_bins) * bin_size;

for i_cluster = 1:nClusters
    % Get all spike times for this cluster
    spike_times = all_spike_times(all_spike_clusters == cluster_ids(i_cluster));

    % Mean firing rate across the session
    if session_duration > 0
        mean_fr_session(i_cluster) = length(spike_times) / session_duration;
    else
        mean_fr_session(i_cluster) = 0;
    end

    % Sparseness: fraction of 100ms bins that are empty
    if ~isempty(spike_times) && n_bins > 0
        spike_counts = histcounts(spike_times, bin_edges);
        empty_bins = sum(spike_counts == 0);
        sparse_fraction(i_cluster) = empty_bins / n_bins;
    else
        sparse_fraction(i_cluster) = 1; % No spikes = 100% sparse
    end
end

fprintf('Session-wide metrics calculated. Mean FR range: %.2f - %.2f sp/s\n', ...
    min(mean_fr_session), max(mean_fr_session));
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
%   (b) Mean FR > 3.5 sp/s
%   (c) NOT significantly activated at opposing locations (180 deg apart)
%   (d) Not sparse (<=70% of 100ms bins empty)
%
% PUTATIVE INTERNEURON:
%   - Mean FR > 3.5 sp/s
%   - Task-modulated (some significant effect somewhere)
%   - Fails to meet "classic" criteria
%
% EXCLUDED:
%   - Mean FR <= 3.5 sp/s OR not task-modulated OR sparse

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
group_mean_frs = cell(1, n_groups);
for i_group = 1:n_groups
    group_mean_frs{i_group} = nan(nClusters, nEpochs);
end

% Test each neuron at each location
for i_cluster = 1:nClusters
    for i_group = 1:n_groups
        trial_mask = trial_groups{i_group};

        if sum(trial_mask) < 2
            continue;
        end

        % Extract firing rates for the current neuron and trial group
        neuron_frs_all_trials = squeeze(epoch_frs(i_cluster, :, trial_mask))';
        neuron_frs = neuron_frs_all_trials(~any(isnan(neuron_frs_all_trials), 2), :);

        if size(neuron_frs, 1) < 2
            continue;
        end

        mean_fr_this_group = mean(neuron_frs, 1, 'omitnan');
        group_mean_frs{i_group}(i_cluster, :) = mean_fr_this_group;

        try
            [p_friedman, ~] = friedman(neuron_frs, 1, 'off');
            if p_friedman < 0.05
                alpha_corr = 0.05 / 3;
                comparisons = [1 2; 1 3; 1 4]; % Bsl vs Vis, Del, Sac

                for i_comp = 1:size(comparisons, 1)
                    p = ranksum(neuron_frs(:, comparisons(i_comp, 1)), ...
                        neuron_frs(:, comparisons(i_comp, 2)));
                    % Significant INCREASE vs baseline
                    if p < alpha_corr && mean(neuron_frs(:, comparisons(i_comp, 2))) > ...
                            mean(neuron_frs(:, comparisons(i_comp, 1)))
                        sig_by_location(i_cluster, i_group, i_comp) = true;
                    end
                end
            end
        catch ME
            fprintf('Stat test failed for cluster %d, group %d: %s\n', ...
                cluster_ids(i_cluster), i_group, ME.message);
        end
    end
end

% Now classify each neuron
for i_cluster = 1:nClusters
    % Check basic criteria
    is_sparse = sparse_fraction(i_cluster) > MAX_SPARSE_FRACTION;
    has_min_fr = mean_fr_session(i_cluster) > MIN_FR_THRESHOLD;

    % Determine if neuron is task-modulated at any location
    is_task_modulated = any(sig_by_location(i_cluster, :, :), 'all');

    % Check for significant activation at contralateral locations
    has_contralateral_activation = false;
    activated_locations = []; % track which locations show activation

    for i_group = 1:n_groups
        if any(sig_by_location(i_cluster, i_group, :))
            theta = location_angles(i_group);
            activated_locations(end+1) = theta;
            if is_contralateral(theta)
                has_contralateral_activation = true;
            end
        end
    end

    % Check for opposing location activation (180 deg apart)
    has_opposing_activation = false;
    for i = 1:length(activated_locations)
        for j = (i+1):length(activated_locations)
            angle_diff = abs(activated_locations(i) - activated_locations(j));
            % Check if locations are approximately 180 deg apart
            if abs(angle_diff - 180) < 30 || abs(angle_diff - 180 + 360) < 30
                has_opposing_activation = true;
                break;
            end
        end
        if has_opposing_activation; break; end
    end

    % Apply classification logic
    if ~has_min_fr || is_sparse || ~is_task_modulated
        % EXCLUDED: low FR, sparse, or not task-modulated
        neuron_class{i_cluster} = 'excluded';
        selected_neurons(i_cluster) = false;
    elseif has_contralateral_activation && ~has_opposing_activation && ~is_sparse
        % CLASSIC: all criteria met
        neuron_class{i_cluster} = 'classic';
        selected_neurons(i_cluster) = true;
    else
        % INTERNEURON: FR > threshold, task-modulated, but fails classic criteria
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
    figFileName = fullfile(output_dir, [session_data.metadata.unique_id, ...
        '_sc_epoch_frs.pdf']);
    pdfSave(figFileName, fig.Position(3:4)/72, fig);
    close(fig);
end

fprintf('Finished screening. Selected %d neurons (Classic: %d, Interneuron: %d).\n', ...
    nnz(selected_neurons), n_classic, n_interneuron);

end
