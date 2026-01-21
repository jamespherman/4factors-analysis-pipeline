function modulation_info = compute_task_modulation(session_data, core_data, condition_defs)
% COMPUTE_TASK_MODULATION Determine if neurons show task-related modulation.
%
%   modulation_info = COMPUTE_TASK_MODULATION(session_data, core_data, condition_defs)
%
%   This function checks each neuron for significant deviation from baseline
%   at ANY target location for ANY task event. Modulation is detected if N
%   consecutive bins (see condition_defs.neuron_inclusion.consecutive_bins_required)
%   show significant deviation (increase OR decrease) from the baseline window.
%
%   Statistical test: arrayROC with 500 bootstrap replicates. Significance is
%   determined by whether the 95% CI excludes 0.5 (chance level). Baseline
%   firing rates are pooled across all bins in the baseline window.
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
%
%   Author: Claude Code
%   Date: 2026-01-18

fprintf('compute_task_modulation: Checking for task modulation...\n');

%% Extract configuration
modulation_windows = condition_defs.modulation_windows;
consecutive_bins_required = condition_defs.neuron_inclusion.consecutive_bins_required;
% Note: Significance is determined by arrayROC bootstrap CI (alpha=0.05 internally)

%% Get neuron and trial information
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nNeurons = numel(cluster_ids);

% Get trial conditions to filter by location
[conditions, ~] = define_task_conditions(session_data);
targetLocIdx = conditions.targetLocIdx;  % Location index (1-4) for each trial

%% Events to check
events_to_check = fieldnames(modulation_windows);

%% Initialize output structure
modulation_info = struct();
for i = 1:nNeurons
    modulation_info(i).cluster_id = cluster_ids(i);
    modulation_info(i).is_modulated = false;
    modulation_info(i).modulated_at = {};
    modulation_info(i).modulation_details = struct();
    modulation_info(i).bin_significance = struct();  % For plotting markers
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
                direction_symbol = '';
                if strcmp(modulation_direction, 'increase')
                    direction_symbol = char(8593);  % Up arrow
                elseif strcmp(modulation_direction, 'decrease')
                    direction_symbol = char(8595);  % Down arrow
                end
                modulated_at_list{end+1} = sprintf('%s@L%d%s', ...
                    event_name, loc, direction_symbol);
            end
        end  % location loop
    end  % event loop

    % Store results for this neuron
    modulation_info(i_neuron).is_modulated = is_modulated_any;
    modulation_info(i_neuron).modulated_at = modulated_at_list;
    modulation_info(i_neuron).modulation_details = details;
    modulation_info(i_neuron).bin_significance = bin_sig;

end  % neuron loop

%% Summary statistics
n_modulated = sum([modulation_info.is_modulated]);
fprintf('compute_task_modulation: %d/%d neurons (%.1f%%) show task modulation.\n', ...
    n_modulated, nNeurons, 100 * n_modulated / nNeurons);

end
