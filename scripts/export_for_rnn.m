function export_for_rnn(session_name, output_path, opts)
% EXPORT_FOR_RNN Export session data for RNN fitting pipeline.
%
% This function loads a processed session_data.mat file and exports it
% in the format specified in docs/DATA_SPEC.md for the Python RNN
% fitting pipeline.
%
% USAGE:
%   export_for_rnn('Feynman_08_15_2025_SC', output_path);
%   export_for_rnn(session_name, output_path, opts);
%
% INPUTS:
%   session_name - String, session identifier (e.g., 'Feynman_08_15_2025_SC')
%   output_path  - String, directory where to save the exported file
%   opts         - (Optional) struct with fields:
%                    .bin_size_ms     - Bin size in milliseconds (default: 25)
%                    .pre_fix_ms      - Pre-fixation buffer (default: 200)
%                    .post_reward_ms  - Post-reward buffer (default: 500)
%                    .data_root       - Path to processed data (default: auto)
%
% OUTPUT:
%   Saves rnn_export_<session_name>.mat to output_path
%
% See also: screen_sc_neurons, screen_snc_neurons, alignAndBinSpikes

%% Input validation and defaults
if nargin < 3
    opts = struct();
end

% Determine project root and add paths
script_path = fileparts(mfilename('fullpath'));
project_root = fileparts(script_path);
addpath(fullfile(project_root, 'code'));
addpath(fullfile(project_root, 'code', 'utils'));

% Default parameters
if ~isfield(opts, 'bin_size_ms')
    opts.bin_size_ms = 25;
end
if ~isfield(opts, 'pre_fix_ms')
    opts.pre_fix_ms = 200;
end
if ~isfield(opts, 'post_reward_ms')
    opts.post_reward_ms = 500;
end

% Default data locations
if ~isfield(opts, 'data_root')
    % Try local cache first, then OneDrive
    local_path = fullfile(project_root, 'data', 'processed', session_name);
    onedrive_path = fullfile(getenv('OneDrive'), 'Neuronal Data Analysis', session_name);

    if exist(local_path, 'dir')
        opts.data_root = local_path;
    elseif exist(onedrive_path, 'dir')
        opts.data_root = onedrive_path;
    else
        error('Cannot find data directory for session: %s', session_name);
    end
end

fprintf('=== RNN Export for session: %s ===\n', session_name);
fprintf('Data root: %s\n', opts.data_root);
fprintf('Bin size: %d ms, Pre-fix: %d ms, Post-reward: %d ms\n', ...
    opts.bin_size_ms, opts.pre_fix_ms, opts.post_reward_ms);

%% 1. Load session data
session_data_path = fullfile(opts.data_root, [session_name '_session_data.mat']);
if ~exist(session_data_path, 'file')
    % Try alternate naming
    session_data_path = fullfile(opts.data_root, 'session_data.mat');
end

if ~exist(session_data_path, 'file')
    error('Session data not found: %s', session_data_path);
end

fprintf('Loading session data from: %s\n', session_data_path);
load(session_data_path, 'session_data');

%% 2. Filter trials to gSac_4factors task only
codes = initCodes();
trialInfo = session_data.trialInfo;
eventTimes = session_data.eventTimes;

% Find valid gSac_4factors trials with successful completion: the trial
% ended in the sacComplete state (21) and reward was delivered.
isGSac4factors = trialInfo.taskCode == codes.uniqueTaskCode_gSac_4factors;
isSuccessful = trialInfo.trialEndState(:) == 21 & ...
    ~cellfun(@isempty, eventTimes.rewardCell(:));
valid_trials = find(isGSac4factors(:) & isSuccessful);

n_trials = length(valid_trials);
fprintf('Found %d valid gSac_4factors trials\n', n_trials);

if n_trials == 0
    error('No valid trials found for session %s', session_name);
end

%% 3. Get neuron classification
[neuron_ids, neuron_type, brain_area_code] = classify_neurons_for_export(session_data, project_root, session_name);

n_neurons = length(neuron_ids);
fprintf('Selected %d neurons (Type 1: %d, Type 2: %d)\n', ...
    n_neurons, sum(neuron_type == 1), sum(neuron_type == 2));

if n_neurons == 0
    error('No neurons passed classification criteria');
end

%% 4. Compute trial timing and bin structure
% All trials are aligned to fixOn
% Time window: -pre_fix_ms to (max_trial_duration + post_reward_ms)

% Convert to seconds
bin_size_s = opts.bin_size_ms / 1000;
pre_fix_s = opts.pre_fix_ms / 1000;
post_reward_s = opts.post_reward_ms / 1000;

% Compute trial durations (fixOn to first reward)
trial_durations = zeros(n_trials, 1);
for t = 1:n_trials
    trial_idx = valid_trials(t);
    fix_on = eventTimes.fixOn(trial_idx);
    reward_times = eventTimes.rewardCell{trial_idx};
    if ~isempty(reward_times)
        trial_durations(t) = reward_times(1) - fix_on;
    else
        % Fallback to trial end if no reward times
        trial_durations(t) = eventTimes.trialEnd(trial_idx) - fix_on;
    end
end

% Sanity check: filter out trials with unreasonable durations
max_reasonable_duration_s = 5.0;  % 5 seconds max per trial
bad_trials = trial_durations > max_reasonable_duration_s | trial_durations <= 0;
if any(bad_trials)
    fprintf('WARNING: Excluding %d trials with invalid durations (<=0 or >%.1fs)\n', ...
        sum(bad_trials), max_reasonable_duration_s);
    fprintf('  Bad durations: min=%.3f, max=%.3f seconds\n', ...
        min(trial_durations(bad_trials)), max(trial_durations(bad_trials)));
    valid_trials = valid_trials(~bad_trials);
    trial_durations = trial_durations(~bad_trials);
    n_trials = length(valid_trials);
    fprintf('  Remaining trials: %d\n', n_trials);
end

max_trial_duration = max(trial_durations);
total_window_s = pre_fix_s + max_trial_duration + post_reward_s;

% Create bin edges
n_time_bins = floor(total_window_s / bin_size_s);
time_axis_s = (-pre_fix_s) + (0:(n_time_bins-1)) * bin_size_s + bin_size_s/2;
time_axis_ms = time_axis_s * 1000;

fprintf('Time window: %.0f to %.0f ms (%d bins)\n', ...
    time_axis_ms(1), time_axis_ms(end), n_time_bins);

%% 5. Bin spike data into firing rates
fprintf('Binning spikes for %d neurons x %d trials...\n', n_neurons, n_trials);

firing_rates = nan(n_neurons, n_time_bins, n_trials);

all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

for n = 1:n_neurons
    cluster_id = neuron_ids(n);
    neuron_spikes = all_spike_times(all_spike_clusters == cluster_id);

    for t = 1:n_trials
        trial_idx = valid_trials(t);
        fix_on = eventTimes.fixOn(trial_idx);

        % Define window for this trial
        win_start = fix_on - pre_fix_s;
        win_end = fix_on + max_trial_duration + post_reward_s;

        % Get spikes in window
        trial_spikes = neuron_spikes(neuron_spikes >= win_start & neuron_spikes < win_end);

        % Align to fixOn
        aligned_spikes = trial_spikes - fix_on;

        % Bin spikes
        bin_edges = (-pre_fix_s) + (0:n_time_bins) * bin_size_s;
        spike_counts = histcounts(aligned_spikes, bin_edges);

        % Convert to firing rate (sp/s)
        firing_rates(n, :, t) = spike_counts / bin_size_s;
    end
end

%% 6. Construct input signals
fprintf('Constructing input signals...\n');

% Initialize all inputs
input_fixation_on = zeros(n_time_bins, n_trials);
input_target_loc = zeros(4, n_time_bins, n_trials);
input_go_signal = zeros(n_time_bins, n_trials);
input_reward_on = zeros(n_time_bins, n_trials);
input_is_face = zeros(n_time_bins, n_trials);
input_is_nonface = zeros(n_time_bins, n_trials);
input_is_bullseye = zeros(n_time_bins, n_trials);
input_high_salience = zeros(n_time_bins, n_trials);
input_low_salience = zeros(n_time_bins, n_trials);

% Eye position - extract from trialInfo.eyeX and eyeY
% These are stored in degrees of visual angle, sampled at display refresh rate (~100 Hz)
input_eye_x = zeros(n_time_bins, n_trials);
input_eye_y = zeros(n_time_bins, n_trials);

% Check if eye position data fields exist
has_eye_position = isfield(trialInfo, 'eyeX') && isfield(trialInfo, 'eyeY') && isfield(trialInfo, 'eyeT');
if ~has_eye_position
    warning('Eye position data (eyeX, eyeY, eyeT) not found in trialInfo. Eye position will be zeros.');
end

for t = 1:n_trials
    trial_idx = valid_trials(t);
    fix_on = eventTimes.fixOn(trial_idx);

    % Get event times relative to fixOn
    target_on_rel = eventTimes.targetOn(trial_idx) - fix_on;
    fix_off_rel = eventTimes.fixOff(trial_idx) - fix_on;
    reward_times = eventTimes.rewardCell{trial_idx};

    % Get trial info
    target_loc = trialInfo.targetLocIdx(trial_idx);
    stim_type = trialInfo.stimType(trial_idx);
    salience = trialInfo.salience(trial_idx);

    % For each time bin, determine signal values
    for b = 1:n_time_bins
        t_bin = time_axis_s(b);  % Time of this bin center relative to fixOn

        % Fixation on: from fixOn until fixOff
        if t_bin >= 0 && t_bin < fix_off_rel
            input_fixation_on(b, t) = 1;
        end

        % Target location: one-hot from targetOn onwards
        if t_bin >= target_on_rel
            input_target_loc(target_loc, b, t) = 1;
        end

        % Go signal: from fixOff onwards
        if t_bin >= fix_off_rel
            input_go_signal(b, t) = 1;
        end

        % Reward: during juice delivery periods
        if ~isempty(reward_times)
            for r_idx = 1:length(reward_times)
                reward_rel = reward_times(r_idx) - fix_on;
                % Assume reward pulse lasts ~200ms
                if t_bin >= reward_rel && t_bin < (reward_rel + 0.2)
                    input_reward_on(b, t) = 1;
                end
            end
        end

        % Target features (revealed at target onset)
        if t_bin >= target_on_rel
            if stim_type == 1
                input_is_face(b, t) = 1;
            elseif stim_type == 2
                input_is_nonface(b, t) = 1;
            elseif stim_type > 2
                input_is_bullseye(b, t) = 1;
                % Salience only applies to bullseye trials
                if salience == 1
                    input_high_salience(b, t) = 1;
                elseif salience == 2
                    input_low_salience(b, t) = 1;
                end
            end
        end
    end

    % Eye position interpolation (if available)
    % Eye data is stored in trialInfo.eyeX, eyeY, eyeT (degrees of visual angle)
    % Eye time is on the VIEWPixx clock - must align to neural clock
    if has_eye_position && ~isempty(trialInfo.eyeX{trial_idx}) && ...
            ~isempty(trialInfo.eyeY{trial_idx}) && ~isempty(trialInfo.eyeT{trial_idx})

        eye_x_raw = trialInfo.eyeX{trial_idx};
        eye_y_raw = trialInfo.eyeY{trial_idx};
        eye_time_raw = trialInfo.eyeT{trial_idx};

        % Align eye time from VIEWPixx clock to neural clock
        % Formula: eye_time_neural = eyeT - pdsTrialStartDP + trialBegin
        pds_trial_start = eventTimes.pdsTrialStartDP(trial_idx);
        trial_begin = eventTimes.trialBegin(trial_idx);
        eye_time_neural = eye_time_raw - pds_trial_start + trial_begin;

        % Convert to time relative to fixOn
        eye_time_rel = eye_time_neural - fix_on;

        % Find valid (non-NaN) eye samples
        valid_eye = ~isnan(eye_x_raw) & ~isnan(eye_y_raw) & ~isnan(eye_time_rel);

        if sum(valid_eye) > 2
            % Interpolate to bin centers
            input_eye_x(:, t) = interp1(eye_time_rel(valid_eye), ...
                eye_x_raw(valid_eye), time_axis_s, 'linear', 0);
            input_eye_y(:, t) = interp1(eye_time_rel(valid_eye), ...
                eye_y_raw(valid_eye), time_axis_s, 'linear', 0);
        end
    end
end

%% 7. Extract trial-level labels
fprintf('Extracting trial-level labels...\n');

trial_reward = zeros(n_trials, 1);
trial_probability = zeros(n_trials, 1);
trial_identity = zeros(n_trials, 1);
trial_salience = zeros(n_trials, 1);
trial_location = zeros(n_trials, 1);
trial_duration_ms = trial_durations * 1000;

% Trial outcome (trialEndState; 21 = success) and the trial's index in
% session_data (also session_trial in data/trial_audit), for tracing.
trial_outcome = trialInfo.trialEndState(valid_trials);
trial_outcome = trial_outcome(:);
trial_session_index = valid_trials(:);

% Each neuron's stable period (from screening's stability detection), as
% session trial indices and as a [n_neurons x n_trials] mask over the
% exported trials. Analyses use only a neuron's stable trials. Neurons
% without stability info are marked stable on all trials.
neuron_stable_range = nan(n_neurons, 2);
neuron_stable_trials = true(n_neurons, n_trials);
if isfield(session_data, 'stability_info')
    stab_ids = [session_data.stability_info.cluster_id];
    for n = 1:n_neurons
        k = find(stab_ids == neuron_ids(n), 1);
        if isempty(k)
            continue;
        end
        si_n = session_data.stability_info(k);
        if si_n.is_stable
            neuron_stable_range(n, :) = [si_n.stable_trial_start, ...
                si_n.stable_trial_end];
            neuron_stable_trials(n, :) = ...
                trial_session_index' >= si_n.stable_trial_start & ...
                trial_session_index' <= si_n.stable_trial_end;
        else
            neuron_stable_trials(n, :) = false;
        end
    end
end

for t = 1:n_trials
    trial_idx = valid_trials(t);

    % Reward: based on rewardDuration threshold (>200ms = high)
    trial_reward(t) = double(trialInfo.rewardDuration(trial_idx) > 200);

    % Probability: block-dependent (Block 1->Loc 1 high, Block 2->Loc 3 high)
    block = trialInfo.blockNumber(trial_idx);
    loc = trialInfo.targetLocIdx(trial_idx);
    if (block == 1 && loc == 1) || (block == 2 && loc == 3)
        trial_probability(t) = 1;  % High probability
    else
        trial_probability(t) = 0;  % Low probability
    end

    % Identity: 1=face, 2=nonface, 3=bullseye
    stim_type = trialInfo.stimType(trial_idx);
    if stim_type == 1
        trial_identity(t) = 1;  % Face
    elseif stim_type == 2
        trial_identity(t) = 2;  % Nonface
    else
        trial_identity(t) = 3;  % Bullseye
    end

    % Salience: 0=N/A (face/nonface), 1=high, 2=low
    salience = trialInfo.salience(trial_idx);
    if stim_type <= 2
        trial_salience(t) = 0;  % N/A for image trials
    else
        trial_salience(t) = salience;  % 1=high, 2=low for bullseye
    end

    % Location: 1-4
    trial_location(t) = loc;
end

%% 8. Prepare metadata
session_name_out = session_name;
export_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
pipeline_version = '1.2.0';

%% 9. Save to output file
if ~exist(output_path, 'dir')
    mkdir(output_path);
end

output_file = fullfile(output_path, sprintf('rnn_export_%s.mat', session_name));
fprintf('Saving to: %s\n', output_file);

% Prepare variables with correct names for export
time_axis = time_axis_ms;          % Time axis in milliseconds
brain_area = brain_area_code;      % Brain area code per neuron
bin_size_ms = opts.bin_size_ms;    % Bin size in ms (scalar)
session_name = session_name_out;   % Session name string

save(output_file, ...
    'firing_rates', 'neuron_ids', 'neuron_type', 'brain_area', ...
    'n_trials', 'n_neurons', 'n_time_bins', ...
    'bin_size_ms', 'time_axis', 'trial_duration_ms', ...
    'input_fixation_on', 'input_target_loc', 'input_go_signal', 'input_reward_on', ...
    'input_eye_x', 'input_eye_y', ...
    'input_is_face', 'input_is_nonface', 'input_is_bullseye', ...
    'input_high_salience', 'input_low_salience', ...
    'trial_reward', 'trial_probability', 'trial_identity', ...
    'trial_salience', 'trial_location', ...
    'trial_outcome', 'trial_session_index', ...
    'neuron_stable_trials', 'neuron_stable_range', ...
    'session_name', 'export_date', 'pipeline_version', ...
    '-v7.3');

fprintf('\n=== Export complete ===\n');
fprintf('Neurons: %d (E=%d, I=%d)\n', n_neurons, sum(neuron_type==1), sum(neuron_type==2));
fprintf('Trials: %d\n', n_trials);
fprintf('Time bins: %d (%.1f to %.1f ms)\n', n_time_bins, time_axis(1), time_axis(end));
fprintf('Output: %s\n', output_file);

end


function [neuron_ids, neuron_type, brain_area_code] = classify_neurons_for_export(session_data, project_root, session_name)
% CLASSIFY_NEURONS_FOR_EXPORT Apply neuron classification for RNN export.
%
% Maps classification results to:
%   Type 1 = Classic/Excitatory (SC 'classic' or SNc 'modulated')
%   Type 2 = Inhibitory (SC 'interneuron')
%
% Excludes neurons that don't pass screening criteria.

% Get brain area from metadata if available, otherwise extract from session name
if isfield(session_data, 'metadata') && isfield(session_data.metadata, 'brain_area')
    brain_area = session_data.metadata.brain_area;
else
    % Extract brain area from session name (e.g., 'Newton_08_14_2025_SC' -> 'SC')
    parts = strsplit(session_name, '_');
    brain_area = parts{end};
    fprintf('Note: No metadata found. Extracted brain_area=''%s'' from session name.\n', brain_area);
end

% Determine brain area code
if strcmpi(brain_area, 'SC')
    brain_area_code_val = 1;
elseif strcmpi(brain_area, 'SNc')
    brain_area_code_val = 2;
else
    brain_area_code_val = 0;
end

% Check if screening has already been run
if isfield(session_data, 'analysis') && isfield(session_data.analysis, 'screening_info')
    % Use existing classification
    screening_info = session_data.analysis.screening_info;
    n_neurons = length(screening_info);

    neuron_ids = [];
    neuron_type = [];

    for i = 1:n_neurons
        class = screening_info(i).neuron_class;

        if strcmpi(brain_area, 'SC')
            if strcmp(class, 'classic')
                neuron_ids(end+1) = screening_info(i).cluster_id;
                neuron_type(end+1) = 1;  % Excitatory
            elseif strcmp(class, 'interneuron')
                neuron_ids(end+1) = screening_info(i).cluster_id;
                neuron_type(end+1) = 2;  % Inhibitory
            end
            % 'excluded' neurons are not included

        elseif strcmpi(brain_area, 'SNc')
            if strcmp(class, 'modulated')
                neuron_ids(end+1) = screening_info(i).cluster_id;
                neuron_type(end+1) = 1;  % All SNc modulated neurons -> excitatory
            end
            % 'excluded' neurons are not included
        end
    end

    neuron_ids = neuron_ids(:);
    neuron_type = neuron_type(:);
    brain_area_code = repmat(brain_area_code_val, length(neuron_ids), 1);

else
    % Need to run screening
    fprintf('No existing screening results found. Running neuron classification...\n');

    % Get condition_defs for screening parameters
    [~, condition_defs] = define_task_conditions();

    cluster_ids = session_data.spikes.cluster_info.cluster_id;
    n_clusters = length(cluster_ids);

    if strcmpi(brain_area, 'SC')
        % Run SC screening
        [selected_neurons, ~, ~, neuron_class] = screen_sc_neurons(session_data, project_root, condition_defs, session_name);

        neuron_ids = [];
        neuron_type = [];

        for i = 1:n_clusters
            if selected_neurons(i)
                neuron_ids(end+1) = cluster_ids(i);
                if strcmp(neuron_class{i}, 'classic')
                    neuron_type(end+1) = 1;
                else  % 'interneuron'
                    neuron_type(end+1) = 2;
                end
            end
        end

    elseif strcmpi(brain_area, 'SNc')
        % For SNc, we need core_data which requires prepare_core_data
        % Simplified approach: use basic FR and sparsity criteria
        fprintf('Running simplified SNc classification (FR > 1 sp/s, sparsity <= 70%%)\n');

        all_spike_times = session_data.spikes.times;
        all_spike_clusters = session_data.spikes.clusters;

        neuron_ids = [];
        neuron_type = [];

        for i = 1:n_clusters
            cluster_id = cluster_ids(i);
            spike_times = all_spike_times(all_spike_clusters == cluster_id);

            if length(spike_times) < 2
                continue;
            end

            % Compute mean FR over active period
            active_duration = spike_times(end) - spike_times(1);
            if active_duration < 1
                mean_fr = length(spike_times);
            else
                mean_fr = length(spike_times) / active_duration;
            end

            % Compute sparsity
            if active_duration >= 1
                edges = spike_times(1):0.1:spike_times(end);
                if length(edges) >= 2
                    counts = histcounts(spike_times, edges);
                    sparsity = sum(counts == 0) / length(counts);
                else
                    sparsity = 0;
                end
            else
                sparsity = 0;
            end

            % Apply thresholds
            if mean_fr >= 1 && sparsity <= 0.7
                neuron_ids(end+1) = cluster_id;
                neuron_type(end+1) = 1;  % All SNc -> excitatory
            end
        end
    else
        error('Unknown brain area: %s', brain_area);
    end

    neuron_ids = neuron_ids(:);
    neuron_type = neuron_type(:);
    brain_area_code = repmat(brain_area_code_val, length(neuron_ids), 1);
end

fprintf('Classification complete: %d neurons selected\n', length(neuron_ids));

end
