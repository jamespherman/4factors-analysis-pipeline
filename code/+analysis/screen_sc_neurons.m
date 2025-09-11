function [selected_neurons, sig_epoch_comparison, scSide, trial_groups, group_sig_results] = screen_sc_neurons(session_data)
% screen_sc_neurons - Implements an inclusive, multi-group method to identify
% task-modulated SC neurons.
%
% This function refactors the neuron selection process to be more inclusive
% and scientifically robust. It calculates firing rates for all memory-guided
% saccade trials from both 'gSac_jph' and 'gSac_4factors' tasks at once.
% It then determines the SC's recorded side hierarchically, preferring
% 'gSac_jph' data. Finally, it tests each neuron for significant modulation
% across multiple, distinct trial groups (all 'gSac_jph' trials, and
% 'gSac_4factors' trials for each unique target location). A neuron is
% selected if it shows significant modulation in ANY of these groups,
% ensuring no task-relevant neurons are missed.
%
% INPUTS:
%   session_data - A struct containing session-specific data, conforming to the
%                  `session_data_dictionary.md`.
%
% OUTPUT:
%   selected_neurons - A logical vector (nClusters x 1) where true indicates
%                      a neuron that passed the selection criteria.
%   sig_epoch_comparison - A logical matrix (nClusters x 3) indicating
%                          significant firing rate changes between epochs.
%   scSide           - A string ('right' or 'left') indicating the determined
%                      recorded SC side.
%   trial_groups     - A cell array of logical masks, each defining a
%                      group of trials used for analysis.
%   group_sig_results - A cell array containing the significance test
%                       results for each neuron within each trial group.
%

fprintf('screen_sc_neurons: Identifying task-modulated neurons...\n');

% Define output directory and filename
project_root = fullfile(findOneDrive, 'Code', ...
    'tokens-analysis-pipeline');
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

%% 4. Hierarchical scSide Determination
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
    left_trials_mask = is_4factors_trial & (thetas_all > 90 & thetas_all < 270)';
    right_trials_mask = is_4factors_trial & (thetas_all < 90 | thetas_all > 270)';

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

%% 5. Inclusive, Multi-Group Neuron Selection
% Iterate through each neuron and test for significant modulation in any of
% several distinct trial groups. A neuron is selected if it passes the
% criteria for any single group.

% Define the trial groups for statistical testing.
% Group 1: All valid gSac_jph memory-guided saccade trials.
% Groups 2-N: gSac_4factors trials, grouped by each unique target location.

% Find logical indices for each task within the combined trial array
is_jph_trial = ismember(all_memSac_trials, gSac_jph_memSac_trials);
is_4factors_trial = ismember(all_memSac_trials, gSac_4factors_memSac_trials);

trial_groups = {};
if any(is_jph_trial)
    trial_groups{end+1} = is_jph_trial;
end

% Get the unique target locations for the 4factors task
unique_locations_4factors = unique(trialInfo.targetTheta( ...
    gSac_4factors_memSac_trials));

% Get target thetas for all combined trials
thetas_all = trialInfo.targetTheta(all_memSac_trials);

for i_loc = 1:length(unique_locations_4factors)
    loc = unique_locations_4factors(i_loc);
    % Create a mask for trials at this location, ONLY for 4factors trials
    loc_mask = (thetas_all == loc) & is_4factors_trial;
    if any(loc_mask)
        trial_groups{end+1} = loc_mask;
    end
end

% Initialize data stores for summary figure
n_groups = length(trial_groups);
group_mean_frs = cell(1, n_groups);
group_sig_results = cell(1, n_groups);
for i_group = 1:n_groups
    group_mean_frs{i_group} = nan(nClusters, nEpochs);
    group_sig_results{i_group} = false(nClusters, 3);
end

% Main loop: iterate through each neuron
for i_cluster = 1:nClusters
    % Nested loop: iterate through each trial group
    for i_group = 1:length(trial_groups)

        trial_mask = trial_groups{i_group};


        % Ensure there are enough trials in the group for statistical tests
        if sum(trial_mask) < 2
            continue;
        end

        % Extract firing rates for the current neuron and trial group
        neuron_frs_all_trials = squeeze(epoch_frs(i_cluster, :, trial_mask))';

        % Remove trials with any NaN epochs
        neuron_frs = neuron_frs_all_trials(~any(isnan(neuron_frs_all_trials), 2), :);

        if size(neuron_frs, 1) < 2
            continue; % Not enough valid trials for this neuron in this group
        end

        is_significant_in_group = false(1, 3);
        mean_fr_this_group = mean(neuron_frs, 1, 'omitnan');

        try
            [p_friedman, ~] = friedman(neuron_frs, 1, 'off');
            if p_friedman < 0.05
                alpha_corr = 0.05 / 3; % Bonferroni for 3 comparisons
                comparisons = [1 2; 1 3; 1 4]; % Bsl vs Vis, Del, Sac

                for i_comp = 1:size(comparisons, 1)
                    p = ranksum(neuron_frs(:, comparisons(i_comp, 1)), ...
                        neuron_frs(:, comparisons(i_comp, 2)));
                    if p < alpha_corr && mean(neuron_frs(:, ...
                            comparisons(i_comp, 2))) > mean(...
                            neuron_frs(:, comparisons(i_comp, 1)))
                        is_significant_in_group(i_comp) = true;
                    end
                end
            end
        catch ME
            fprintf('Stat test failed for cluster %d, group %d: %s\n', ...
                cluster_ids(i_cluster), i_group, ME.message);
        end

        % Store the results for this group for later plotting
        group_mean_frs{i_group}(i_cluster, :) = mean_fr_this_group;
        group_sig_results{i_group}(i_cluster, :) = is_significant_in_group;

        % If significant, update the master selection and sig results
        if any(is_significant_in_group) && max(mean_fr_this_group) > 5
            selected_neurons(i_cluster) = true;

            % Store the significance profile from the first group that
            % passes the test as the canonical result for the neuron.
            if ~any(sig_epoch_comparison(i_cluster, :))
                sig_epoch_comparison(i_cluster, :) = ...
                    is_significant_in_group;
            end
        end
    end % end of group loop
end % end of cluster loop

% --- Generate new summary figure ---
if n_groups > 0
    % Calculate global max firing rate for consistent color scaling
    global_max_fr = 0;
    for i_group = 1:n_groups
        max_in_group = max(group_mean_frs{i_group}, [], 'all', 'omitnan');
        if max_in_group > global_max_fr
            global_max_fr = max_in_group;
        end
    end
    if global_max_fr == 0; global_max_fr = 1; end % Avoid Clim = [0 0]

    fig = figure('Color', 'w', 'Position', [100, 100, 350 * n_groups, 700]);

    % Define labels for plots
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

        % Generate title for the pair of plots
        title_str = sprintf('Group %d', i_group);
        if i_group == 1 && any(is_jph_trial)
            title_str = 'gSac_jph';
        else
            % Find which location this is for 4factors
            group_mask = trial_groups{i_group};
            theta_for_group = unique(thetas_all(group_mask));
            if ~isempty(theta_for_group)
                title_str = sprintf('gSac_4factors: Theta %d', ...
                    theta_for_group(1)/10);
            end
        end
        title(title_str);

        % Subplot for Significance
        plot_idx_sig = (i_group - 1) * 2 + 2;
        mySubPlot([1, n_groups * 2, plot_idx_sig], 'Width', 0.92, ...
            'LeftMargin', 0.05);
        imagesc(group_sig_results{i_group});
        colormap(gca, flipud(bone));
        set(gca, 'XTick', 1:3, 'XTickLabel', sig_x_labels, ...
            'YTickLabel', []);
    end

    % Save the figure
    figFileName = fullfile(output_dir, [session_data.metadata.unique_id, ...
        '_sc_epoch_frs.pdf']);
    pdfSave(figFileName, fig.Position(3:4)/72, fig);
end

% close the figure window
close(fig);

fprintf('Finished screening. Found %d task-modulated SC neurons.\n', ...
    nnz(selected_neurons));

end
