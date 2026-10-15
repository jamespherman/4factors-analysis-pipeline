function stability_info = detect_stable_trial_period(spike_times, spike_clusters, ...
    trial_start_times, trial_end_times, cluster_ids, params)
% DETECT_STABLE_TRIAL_PERIOD Identify stable recording period for each neuron.
%
%   stability_info = DETECT_STABLE_TRIAL_PERIOD(spike_times, spike_clusters, ...
%       trial_start_times, trial_end_times, cluster_ids, params)
%
%   This function uses a sliding window approach to detect the longest
%   contiguous period where each neuron's firing rate is stable. Inspired by
%   SpikeInterface's presence ratio and SpikeSift's dynamic segmentation.
%
%   ALGORITHM:
%   1. For each neuron, compute FR in sliding windows of N trials
%   2. Calculate local coefficient of variation (CV) of FR
%   3. Find longest run of consecutive windows with CV < threshold
%   4. Return trial indices and metrics for the stable period
%
%   INPUTS:
%       spike_times      - [nSpikes x 1] spike times in seconds
%       spike_clusters   - [nSpikes x 1] cluster ID for each spike
%       trial_start_times- [nTrials x 1] start time of each trial
%       trial_end_times  - [nTrials x 1] end time of each trial
%       cluster_ids      - [nClusters x 1] cluster IDs to analyze
%       params           - struct with fields:
%           .window_size     - Number of trials per window (default: 20)
%           .step_size       - Window step in trials (default: 5)
%           .cv_threshold    - Max CV for "stable" (default: 0.5)
%           .min_fr_threshold- Minimum FR to consider present (default: 0.1 sp/s)
%           .min_stable_trials- Minimum trials in stable period (default: 50)
%
%   OUTPUTS:
%       stability_info   - [nClusters x 1] struct array with fields:
%           .cluster_id          - Cluster ID
%           .stable_trial_start  - First trial index of stable period
%           .stable_trial_end    - Last trial index of stable period
%           .n_stable_trials     - Number of trials in stable period
%           .stable_fr           - Mean FR over stable period (sp/s)
%           .stable_sparsity     - Sparsity over stable period
%           .presence_ratio      - Fraction of windows with firing
%           .is_stable           - Boolean: meets stability criteria
%           .stability_score     - Continuous score (0-1, higher = more stable)
%           .window_frs          - [nWindows x 1] FR in each window
%           .window_cvs          - [nWindows x 1] CV in each window
%           .exclusion_reason    - If not stable, why
%
%   References:
%       SpikeInterface Quality Metrics: https://spikeinterface.readthedocs.io/
%       SpikeSift (2025): https://arxiv.org/html/2504.01604v1
%
%   Author: Claude Code
%   Date: 2026-01-23

% Default parameters
if nargin < 6 || isempty(params)
    params = struct();
end
params = set_defaults(params);

nClusters = length(cluster_ids);
nTrials = length(trial_start_times);

fprintf('detect_stable_trial_period: Analyzing %d neurons across %d trials...\n', ...
    nClusters, nTrials);
fprintf('  Parameters: window_size=%d, step_size=%d, cv_threshold=%.2f, min_stable_trials=%d\n', ...
    params.window_size, params.step_size, params.cv_threshold, params.min_stable_trials);

% Initialize output struct array
stability_info = struct(...
    'cluster_id', num2cell(cluster_ids(:)), ...
    'stable_trial_start', num2cell(ones(nClusters, 1)), ...
    'stable_trial_end', num2cell(nTrials * ones(nClusters, 1)), ...
    'n_stable_trials', num2cell(nTrials * ones(nClusters, 1)), ...
    'stable_fr', num2cell(zeros(nClusters, 1)), ...
    'stable_sparsity', num2cell(zeros(nClusters, 1)), ...
    'presence_ratio', num2cell(zeros(nClusters, 1)), ...
    'is_stable', num2cell(false(nClusters, 1)), ...
    'stability_score', num2cell(zeros(nClusters, 1)), ...
    'window_frs', cell(nClusters, 1), ...
    'window_cvs', cell(nClusters, 1), ...
    'exclusion_reason', cell(nClusters, 1));

% Handle edge case: not enough trials for windowed analysis
if nTrials < params.window_size
    fprintf('detect_stable_trial_period: Only %d trials (< window_size=%d). Using all trials.\n', ...
        nTrials, params.window_size);
    for i_cluster = 1:nClusters
        stability_info(i_cluster).is_stable = true;
        stability_info(i_cluster).exclusion_reason = 'Too few trials for windowed analysis';
    end
    return;
end

% Define sliding windows
window_starts = 1:params.step_size:(nTrials - params.window_size + 1);
nWindows = length(window_starts);

if nWindows < 3
    fprintf('detect_stable_trial_period: Only %d windows (need >= 3 for CV). Using all trials.\n', nWindows);
    for i_cluster = 1:nClusters
        stability_info(i_cluster).is_stable = true;
        stability_info(i_cluster).exclusion_reason = 'Too few windows for CV analysis';
    end
    return;
end

% Process each cluster
n_stable = 0;
for i_cluster = 1:nClusters
    cluster_id = cluster_ids(i_cluster);
    cluster_spike_times = spike_times(spike_clusters == cluster_id);

    % Progress update
    if mod(i_cluster, 20) == 1 || i_cluster == nClusters
        fprintf('detect_stable_trial_period: Processing cluster %d/%d...\n', ...
            i_cluster, nClusters);
    end

    if isempty(cluster_spike_times)
        stability_info(i_cluster).exclusion_reason = 'No spikes';
        stability_info(i_cluster).is_stable = false;
        stability_info(i_cluster).stable_fr = 0;
        continue;
    end

    % Compute FR in each sliding window
    window_frs = zeros(nWindows, 1);
    for i_win = 1:nWindows
        win_trial_start = window_starts(i_win);
        win_trial_end = window_starts(i_win) + params.window_size - 1;

        win_start_time = trial_start_times(win_trial_start);
        win_end_time = trial_end_times(win_trial_end);
        win_duration = win_end_time - win_start_time;

        win_spikes = sum(cluster_spike_times >= win_start_time & ...
                        cluster_spike_times < win_end_time);
        window_frs(i_win) = win_spikes / win_duration;
    end

    stability_info(i_cluster).window_frs = window_frs;

    % Calculate presence ratio (fraction of windows with meaningful firing)
    presence_ratio = mean(window_frs >= params.min_fr_threshold);
    stability_info(i_cluster).presence_ratio = presence_ratio;

    % Calculate local CV using rolling window of 3 consecutive windows
    window_cvs = zeros(nWindows - 2, 1);
    for i_win = 1:(nWindows - 2)
        local_frs = window_frs(i_win:(i_win + 2));
        if mean(local_frs) > 0
            window_cvs(i_win) = std(local_frs) / mean(local_frs);
        else
            window_cvs(i_win) = inf;  % No firing = unstable
        end
    end
    stability_info(i_cluster).window_cvs = window_cvs;

    % Find longest stable run (consecutive windows with CV < threshold)
    is_stable_window = window_cvs < params.cv_threshold;
    [longest_run_start, longest_run_length] = find_longest_run(is_stable_window);

    if longest_run_length > 0
        % Convert window indices to trial indices
        stable_win_start = longest_run_start;
        stable_win_end = longest_run_start + longest_run_length - 1;

        % The stable window indices refer to the CV array, which is offset by 1
        % from the original window_starts (since CV requires 3 consecutive windows)
        stable_trial_start = window_starts(stable_win_start);
        stable_trial_end = min(window_starts(stable_win_end + 2) + params.window_size - 1, nTrials);
        n_stable_trials = stable_trial_end - stable_trial_start + 1;

        % Check minimum trial requirement
        if n_stable_trials >= params.min_stable_trials
            stability_info(i_cluster).stable_trial_start = stable_trial_start;
            stability_info(i_cluster).stable_trial_end = stable_trial_end;
            stability_info(i_cluster).n_stable_trials = n_stable_trials;
            stability_info(i_cluster).is_stable = true;
            n_stable = n_stable + 1;

            % Compute FR and sparsity over stable period only
            stable_start_time = trial_start_times(stable_trial_start);
            stable_end_time = trial_end_times(stable_trial_end);
            stable_spikes = cluster_spike_times(cluster_spike_times >= stable_start_time & ...
                                                cluster_spike_times < stable_end_time);
            stable_duration = stable_end_time - stable_start_time;

            stability_info(i_cluster).stable_fr = length(stable_spikes) / stable_duration;

            % Sparsity over stable period (fraction of 100ms bins empty)
            bin_edges = stable_start_time:0.1:stable_end_time;
            if length(bin_edges) > 1
                counts = histcounts(stable_spikes, bin_edges);
                stability_info(i_cluster).stable_sparsity = sum(counts == 0) / length(counts);
            end

            % Stability score (0-1): combines presence ratio and run length
            stability_info(i_cluster).stability_score = ...
                presence_ratio * (n_stable_trials / nTrials);
        else
            stability_info(i_cluster).exclusion_reason = sprintf(...
                'Stable period too short (%d < %d trials)', ...
                n_stable_trials, params.min_stable_trials);
            stability_info(i_cluster).is_stable = false;

            % Still set the stable period info (for reference)
            stability_info(i_cluster).stable_trial_start = stable_trial_start;
            stability_info(i_cluster).stable_trial_end = stable_trial_end;
            stability_info(i_cluster).n_stable_trials = n_stable_trials;
        end
    else
        stability_info(i_cluster).exclusion_reason = 'No stable period found (CV always high)';
        stability_info(i_cluster).is_stable = false;
    end

    % If not stable, compute FR over all trials as fallback
    if ~stability_info(i_cluster).is_stable
        total_duration = trial_end_times(end) - trial_start_times(1);
        stability_info(i_cluster).stable_fr = length(cluster_spike_times) / total_duration;
    end
end

fprintf('detect_stable_trial_period: %d/%d neurons (%.1f%%) have stable periods.\n', ...
    n_stable, nClusters, 100 * n_stable / nClusters);

end

%% Helper functions
function params = set_defaults(params)
    if ~isfield(params, 'window_size'), params.window_size = 20; end
    if ~isfield(params, 'step_size'), params.step_size = 5; end
    if ~isfield(params, 'cv_threshold'), params.cv_threshold = 0.5; end
    if ~isfield(params, 'min_fr_threshold'), params.min_fr_threshold = 0.1; end
    if ~isfield(params, 'min_stable_trials'), params.min_stable_trials = 50; end
end

function [start_idx, run_length] = find_longest_run(bool_array)
    % Find the longest consecutive run of true values
    if isempty(bool_array) || ~any(bool_array)
        start_idx = 0;
        run_length = 0;
        return;
    end

    % Find run starts and ends
    bool_array = bool_array(:);
    d = diff([0; bool_array; 0]);
    run_starts = find(d == 1);
    run_ends = find(d == -1) - 1;
    run_lengths = run_ends - run_starts + 1;

    [run_length, max_idx] = max(run_lengths);
    start_idx = run_starts(max_idx);
end
