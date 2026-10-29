function modulation_info = compute_task_modulation(session_data, core_data, condition_defs, stability_info)
% COMPUTE_TASK_MODULATION Determine if neurons show task-related modulation.
%
%   modulation_info = COMPUTE_TASK_MODULATION(session_data, core_data, condition_defs)
%
%   This function checks each neuron for significant deviation from baseline
%   at ANY target location for ANY task event. Modulation is detected if N
%   consecutive bins (see condition_defs.neuron_inclusion.consecutive_bins_required)
%   show significant deviation (increase OR decrease) from the baseline window.
%
%   Statistic: per 25 ms test bin, the ROC area of the test-bin rates
%   (across trials) against the baseline rates pooled across trials and
%   baseline bins. Significance depends on
%   condition_defs.neuron_inclusion.modulation_test:
%
%   'random_time' (default, 2026-10-28): the same run score, tested
%       against a random-time null (shared SC spec, section 4): each
%       trial is moved by a random shift within +-5 min (one shift for
%       all its events), which keeps the trial count, event geometry
%       and the neuron's burst structure but no task events. Each event
%       x location x direction is scaled by its own null 99th
%       percentile, and the largest scaled score sets a family-wise p
%       per neuron (single-step max-T; modulation_random_time_null).
%   'permutation' (2026-10-17): a run of N consecutive bins
%       deviating in the same direction is scored by its smallest
%       deviation |AUC - 0.5|. Each neuron's null distribution of the
%       largest run score over all events x locations comes from
%       permuting bins within trials (modulation_run_stat), and the
%       neuron is modulated if p < modulation_alpha. This holds false
%       positives at alpha at any firing rate (the legacy test flagged
%       ~90% of unmodulated neurons at 0.5 sp/s, where most bins are 0).
%   'bootstrap_ci' (legacy): arrayROC with 500 bootstrap replicates; a
%       bin is significant if its 95% CI excludes 0.5, uncorrected.
%
%   IMPORTANT: Modulation is checked separately for each location (1-4) to
%   avoid washing out spatially selective responses through pooling.
%
%   INPUTS:
%       session_data   - Standard session data structure
%       core_data      - Binned spike data from prepare_core_data
%       condition_defs - Analysis configuration with modulation_windows
%
%   OUTPUTS:
%       modulation_info - Struct array (one per neuron) with fields:
%           .cluster_id      - Cluster ID for this neuron
%           .is_modulated    - Logical: true if modulated at any location/event
%           .modulated_at    - Cell array of strings describing where modulation
%                              was detected, e.g., {'targetOn_loc1', 'saccadeOnset_loc3'}
%           .modulation_details - Struct with detailed location x event results
%           .bin_significance   - Struct with per-bin significance for plotting
%                                 Format: .{event}.{loc}.above, .below, .time_vector
%           .modulation_p    - Neuron-level family-wise p (NaN for
%                              legacy)
%           .method          - modulation_test used (cache key)
%
%   Author: Claude Code
%   Date: 2026-01-18

fprintf('compute_task_modulation: Checking for task modulation...\n');

% Handle optional stability_info argument for neuron present trials
if nargin < 4
    stability_info = [];
end

%% Extract configuration
modulation_windows = condition_defs.modulation_windows;
neuron_inclusion = condition_defs.neuron_inclusion;
consecutive_bins_required = neuron_inclusion.consecutive_bins_required;
method = 'random_time';
if isfield(neuron_inclusion, 'modulation_test')
    method = neuron_inclusion.modulation_test;
end
use_perm = strcmp(method, 'permutation');
use_rt = strcmp(method, 'random_time');
use_blocks = use_perm || use_rt;
n_perm = 1000;
if isfield(neuron_inclusion, 'modulation_n_perm')
    n_perm = neuron_inclusion.modulation_n_perm;
end
alpha = 0.05;
if isfield(neuron_inclusion, 'modulation_alpha')
    alpha = neuron_inclusion.modulation_alpha;
end
fprintf('compute_task_modulation: method = %s\n', method);

% Permutations are seeded per neuron (by cluster ID) so results do not
% depend on processing order; the caller's random state is restored
rng_state = rng;

%% Get neuron and trial information
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nNeurons = numel(cluster_ids);

% Get trial conditions to filter by location
[conditions, ~] = define_task_conditions(session_data);
targetLocIdx = conditions.targetLocIdx;  % Location index (1-4) for each trial

%% Compute neuron present trial masks if stability_info is provided
use_stable_trials = false;
neuron_trial_masks = [];
if ~isempty(stability_info) && isfield(core_data, 'trial_indices')
    neuron_trial_masks = get_neuron_present_trials(stability_info, ...
        core_data.trial_indices, nNeurons);
    use_stable_trials = true;
    fprintf('compute_task_modulation: Using neuron present trials (stable periods)\n');
end

%% Events to check
events_to_check = fieldnames(modulation_windows);

%% Random-time null: event times of every core trial, the windows that
% must stay inside the recording, and the recording span
if use_rt
    rt_opts = struct('bin_width', 0.025, ...
        'n_run', consecutive_bins_required, 'n_null', n_perm, ...
        'max_shift_s', 300, 'scale_floor', 0.02, ...
        'rec_span', [min(session_data.spikes.times), ...
        max(session_data.spikes.times)]);
    if isfield(neuron_inclusion, 'modulation_max_shift_s')
        rt_opts.max_shift_s = neuron_inclusion.modulation_max_shift_s;
    end
    tri = core_data.trial_indices(:);
    trial_events = nan(numel(tri), numel(events_to_check));
    rt_opts.event_windows = zeros(numel(events_to_check), 2);
    for i_event = 1:numel(events_to_check)
        ev_name = events_to_check{i_event};
        if strcmp(ev_name, 'reward')
            % first reward pulse, as in prepare_neuronal_data
            rc = session_data.eventTimes.rewardCell(tri);
            first = nan(numel(tri), 1);
            has = ~cellfun(@isempty, rc);
            first(has) = cellfun(@(c) c(1), rc(has));
            trial_events(:, i_event) = first;
        else
            trial_events(:, i_event) = ...
                session_data.eventTimes.(ev_name)(tri);
        end
        rt_opts.event_windows(i_event, :) = ...
            condition_defs.event_windows.(ev_name);
    end
end

%% Initialize output structure
modulation_info = struct();
for i = 1:nNeurons
    modulation_info(i).cluster_id = cluster_ids(i);
    modulation_info(i).is_modulated = false;
    modulation_info(i).modulated_at = {};
    modulation_info(i).modulation_details = struct();
    modulation_info(i).bin_significance = struct();  % For plotting markers
    modulation_info(i).modulation_p = NaN;
    modulation_info(i).method = method;
end

%% Main loop: check each neuron
for i_neuron = 1:nNeurons

    % Progress update every 10 neurons
    if mod(i_neuron, 10) == 1 || i_neuron == nNeurons
        fprintf('compute_task_modulation: Processing neuron %d/%d...\n', i_neuron, nNeurons);
    end

    is_modulated_any = false;
    modulated_at_list = {};
    details = struct();
    bin_sig = struct();  % Store per-bin significance for plotting

    % Permutation test: run scores of every event x location, and the
    % largest null run score per permutation across all of them
    if use_blocks
        rng(double(cluster_ids(i_neuron)), 'twister');
        blocks = struct('event', {}, 'loc', {}, 'run_obs', {}, ...
            'test_bins', {}, 'event_idx', {}, 'trials', {}, ...
            'bin_starts', {}, 'n_baseline', {});
        null_max = -inf(n_perm, 1);
    end

    % Check each event
    for i_event = 1:length(events_to_check)
        event_name = events_to_check{i_event};

        % Skip if this event is not in core_data
        if ~isfield(core_data.spikes, event_name)
            continue;
        end

        % Get window definitions for this event
        win_def = modulation_windows.(event_name);
        baseline_window = win_def.baseline;  % [start, end] in seconds
        test_window = win_def.test;          % [start, end] in seconds

        % Get spike rate data for this event
        event_data = core_data.spikes.(event_name);
        time_vector = event_data.time_vector;
        rates = event_data.rates;  % [nNeurons x nTrials x nTimeBins]
        bin_width = median(diff(time_vector));

        % Find bin indices for baseline and test windows
        baseline_bins = find(time_vector >= baseline_window(1) & ...
                            time_vector < baseline_window(2));
        test_bins = find(time_vector >= test_window(1) & ...
                        time_vector < test_window(2));

        if isempty(baseline_bins) || isempty(test_bins)
            warning('compute_task_modulation: No bins found for %s. Skipping.', event_name);
            continue;
        end

        % Initialize bin_sig structure for this event
        bin_sig.(event_name) = struct();

        % Check each target location separately
        for loc = 1:4
            loc_trials = find(targetLocIdx == loc);
            loc_key = sprintf('loc%d', loc);

            % Initialize bin significance for this location
            bin_sig.(event_name).(loc_key) = struct(...
                'above', false(1, length(time_vector)), ...
                'below', false(1, length(time_vector)), ...
                'time_vector', time_vector, ...
                'test_window', test_window);

            % Apply neuron-specific stable trial mask (if available)
            if use_stable_trials && ~isempty(neuron_trial_masks)
                stable_mask = neuron_trial_masks(i_neuron, :)';
                loc_trials = loc_trials(stable_mask(loc_trials));
            end

            if length(loc_trials) < 5
                % Not enough trials at this location
                details.(event_name).(loc_key) = struct(...
                    'checked', false, ...
                    'reason', 'insufficient_trials', ...
                    'n_trials', length(loc_trials));
                continue;
            end

            % Extract rates for this neuron at this location
            neuron_rates = squeeze(rates(i_neuron, loc_trials, :));  % [nTrials x nTimeBins]

            % Check for NaN trials and remove them
            valid_trials = ~any(isnan(neuron_rates), 2);
            neuron_rates = neuron_rates(valid_trials, :);

            if size(neuron_rates, 1) < 5
                details.(event_name).(loc_key) = struct(...
                    'checked', false, ...
                    'reason', 'insufficient_valid_trials', ...
                    'n_valid_trials', size(neuron_rates, 1));
                continue;
            end

            if use_blocks
                % Spike counts: baseline bins first, then test bins
                counts = round(neuron_rates(:, ...
                    [baseline_bins(:); test_bins(:)]) * bin_width);
                if use_perm
                    [run_obs, block_null, auc] = modulation_run_stat( ...
                        counts, numel(baseline_bins), ...
                        consecutive_bins_required, n_perm);
                    null_max = max(null_max, block_null);
                else
                    % scored after the loops, with the random-time null
                    run_obs = []; auc = [];
                end
                blocks(end + 1) = struct('event', event_name, ...
                    'loc', loc, 'run_obs', run_obs, ...
                    'test_bins', test_bins(:)', ...
                    'event_idx', i_event, ...
                    'trials', loc_trials(valid_trials), ...
                    'bin_starts', time_vector( ...
                    [baseline_bins(:); test_bins(:)]) - bin_width / 2, ...
                    'n_baseline', numel(baseline_bins)); %#ok<AGROW>
                details.(event_name).(loc_key) = struct(...
                    'checked', true, ...
                    'is_modulated', false, ...
                    'direction', '', ...
                    'auc', auc, ...
                    'p', NaN, ...
                    'n_trials', sum(valid_trials));
                % Significance is assigned after all events and
                % locations, against the neuron's null distribution
                continue;
            end

            % Get baseline firing rate distribution (pooled across baseline bins)
            % Following analyze_baseline_comparison.m approach
            baseline_rates = neuron_rates(:, baseline_bins);  % [nTrials x nBaselineBins]
            n_test_bins = length(test_bins);

            % Flatten baseline and replicate to match test bins
            % This pools all baseline bins together for comparison
            baseline_pooled = repmat(baseline_rates(:), 1, n_test_bins);

            % Get test window rates
            test_rates = neuron_rates(:, test_bins);  % [nTrials x nTestBins]

            % Use arrayROC with bootstrap (500 replicates) for significance
            % sig output: +1 (above chance), -1 (below chance), 0 (not significant)
            % 'notpfp' flag skips ROC curve point computation for speed
            try
                [~, ~, ~, sig_vals] = arrayROC(baseline_pooled, test_rates, 500, 0.05, 'notpfp');
            catch
                sig_vals = zeros(1, n_test_bins);
            end

            % Convert arrayROC sig output to our format
            % sig_vals: +1 means test > baseline, -1 means test < baseline
            direction = sig_vals;  % +1, -1, or 0
            is_sig = sig_vals ~= 0;  % True if significantly different from chance

            % Store bin-by-bin significance for plotting
            % Map test window bins back to full time vector
            sig_above = false(1, length(time_vector));
            sig_below = false(1, length(time_vector));
            for i_bin = 1:n_test_bins
                bin_idx = test_bins(i_bin);
                if is_sig(i_bin)
                    if direction(i_bin) == 1
                        sig_above(bin_idx) = true;
                    else
                        sig_below(bin_idx) = true;
                    end
                end
            end
            bin_sig.(event_name).(loc_key).above = sig_above;
            bin_sig.(event_name).(loc_key).below = sig_below;

            % Find runs of consecutive significant bins
            modulated_this_loc_event = false;
            modulation_direction = '';

            % Check for consecutive significant increases
            for i_start = 1:(n_test_bins - consecutive_bins_required + 1)
                bin_range = i_start:(i_start + consecutive_bins_required - 1);
                if all(is_sig(bin_range)) && all(direction(bin_range) == 1)
                    modulated_this_loc_event = true;
                    modulation_direction = 'increase';
                    break;
                end
                if all(is_sig(bin_range)) && all(direction(bin_range) == -1)
                    modulated_this_loc_event = true;
                    modulation_direction = 'decrease';
                    break;
                end
            end

            % Store details
            details.(event_name).(loc_key) = struct(...
                'checked', true, ...
                'is_modulated', modulated_this_loc_event, ...
                'direction', modulation_direction, ...
                'sig_vals', sig_vals, ...
                'n_trials', sum(valid_trials));

            % Update master flags
            if modulated_this_loc_event
                is_modulated_any = true;
                modulated_at_list{end+1} = format_modulated_at( ...
                    event_name, loc, modulation_direction); %#ok<AGROW>
            end
        end  % location loop
    end  % event loop

    %% Score each event x location against the neuron's null
    if use_blocks && ~isempty(blocks)
        n_run = consecutive_bins_required;
        if use_rt
            cid = double(cluster_ids(i_neuron));
            rt_opts.seed = cid;
            rt = modulation_random_time_null( ...
                session_data.spikes.times( ...
                session_data.spikes.clusters == cid), ...
                blocks, trial_events, rt_opts);
        else
            % p of a run score = fraction of permutations whose
            % largest run score (anywhere) is at least as large
            run_p = @(v) (1 + sum(null_max >= v)) / (1 + n_perm);
        end
        best_p = 1;
        for i_block = 1:numel(blocks)
            blk = blocks(i_block);
            loc_key = sprintf('loc%d', blk.loc);
            if use_rt
                p_block = rt.p_block(i_block);
                best_dir = rt.dir_block(i_block);
                details.(blk.event).(loc_key).auc = rt.auc{i_block};
            else
                % strongest run in each direction (1 = increase)
                [best_score, best_dir] = max(max(blk.run_obs, [], 2));
                p_block = run_p(best_score);
            end
            best_p = min(best_p, p_block);
            details.(blk.event).(loc_key).p = p_block;
            if p_block >= alpha
                continue;
            end
            dir_names = {'increase', 'decrease'};
            details.(blk.event).(loc_key).is_modulated = true;
            details.(blk.event).(loc_key).direction = ...
                dir_names{best_dir};
            is_modulated_any = true;
            modulated_at_list{end+1} = format_modulated_at( ...
                blk.event, blk.loc, dir_names{best_dir}); %#ok<AGROW>
            % Mark the bins of every significant run for plotting
            n_tv = numel(bin_sig.(blk.event).(loc_key).above);
            marks = false(2, n_tv);
            if use_rt
                rp = rt.run_p{i_block};
            else
                rp = arrayfun(run_p, blk.run_obs);
            end
            for i_dir = 1:2
                for i_start = 1:size(rp, 2)
                    if rp(i_dir, i_start) < alpha
                        idx = blk.test_bins(i_start + (0:(n_run - 1)));
                        marks(i_dir, idx) = true;
                    end
                end
            end
            bin_sig.(blk.event).(loc_key).above = marks(1, :);
            bin_sig.(blk.event).(loc_key).below = marks(2, :);
        end
        modulation_info(i_neuron).modulation_p = best_p;
    end

    % Store results for this neuron
    modulation_info(i_neuron).is_modulated = is_modulated_any;
    modulation_info(i_neuron).modulated_at = modulated_at_list;
    modulation_info(i_neuron).modulation_details = details;
    modulation_info(i_neuron).bin_significance = bin_sig;

end  % neuron loop

rng(rng_state);

%% Summary statistics
n_modulated = sum([modulation_info.is_modulated]);
fprintf('compute_task_modulation: %d/%d neurons (%.1f%%) show task modulation.\n', ...
    n_modulated, nNeurons, 100 * n_modulated / nNeurons);

end

function s = format_modulated_at(event_name, loc, direction)
% e.g. 'reward@L2' followed by an up or down arrow
direction_symbol = '';
if strcmp(direction, 'increase')
    direction_symbol = char(8593);  % Up arrow
elseif strcmp(direction, 'decrease')
    direction_symbol = char(8595);  % Down arrow
end
s = sprintf('%s@L%d%s', event_name, loc, direction_symbol);
end
