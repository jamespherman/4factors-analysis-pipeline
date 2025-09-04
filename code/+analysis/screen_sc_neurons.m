function selected_neurons = screen_sc_neurons(session_data, scSide)
% screen_sc_neurons - Identifies task-modulated neurons in the SC.
%
% This function selects neurons that show significant changes in firing rate
% across different task epochs (Baseline, Visual, Delay, Saccade). It uses
% the standardized `session_data` structure.
%
% INPUTS:
%   session_data - A struct containing session-specific data, conforming to the
%                  `session_data_dictionary.md`.
%   scSide       - A string ('right' or 'left') indicating the recorded SC side.
%                  This info typically comes from `config/session_manifest.csv`.
%
% OUTPUT:
%   selected_neurons - A logical vector (nClusters x 1) where true indicates
%                      a neuron that passed the selection criteria.
%

fprintf('screen_sc_neurons: Identifying task-modulated neurons...\n');

% --- Setup and Data Extraction ---
% Get cluster information
cluster_info = session_data.spikes.cluster_info;
cluster_ids = cluster_info.cluster_id;
nClusters = height(cluster_info);

% Get spike times and cluster assignments
all_spike_times = session_data.spikes.times;
all_spike_clusters = session_data.spikes.clusters;

% Get trial information and event times
trialInfo = session_data.trialInfo;
eventTimes = session_data.eventTimes;

% --- Trial Selection ---
% Select high-saliency trials for contralateral targets.
% Assumption: salience == 1 means high saliency.
% Assumption: targTheta > 0 is Left Visual Field, targTheta < 0 is Right.
is_high_saliency = trialInfo.salience == 1;

if strcmpi(scSide, 'right')      % Right SC -> Contralateral is Left VF (targTheta > 0)
    is_contralateral = trialInfo.targTheta > 0;
    fprintf('screen_sc_neurons: Right SC. Contralateral = Left Visual Field (targTheta > 0).\n');
elseif strcmpi(scSide, 'left') % Left SC -> Contralateral is Right VF (targTheta < 0)
    is_contralateral = trialInfo.targTheta < 0;
    fprintf('screen_sc_neurons: Left SC. Contralateral = Right Visual Field (targTheta < 0).\n');
else
    error('Invalid scSide: "%s". Must be "left" or "right".', scSide);
end

trialInd = find(is_high_saliency & is_contralateral);

% Handle edge cases
if isempty(trialInd)
    fprintf('WARNING in screen_sc_neurons: No trials match selection criteria. Skipping.\n');
    selected_neurons = false(nClusters, 1);
    return;
end
if nClusters == 0
    fprintf('WARNING in screen_sc_neurons: No clusters found.\n');
    selected_neurons = false(0, 1);
    return;
end

nTrials = numel(trialInd);
epoch_frs = zeros(nClusters, 4, nTrials); % [nClusters x nEpochs x nTrials]

% --- Epoch Definition ---
% Define analysis epochs relative to key trial events from `session_data.eventTimes`
% {event_field, start_offset, end_offset, duration}
epochs = {
    'targOn',  -0.075, 0.025,  0.1;   % 1. Baseline (rel. to Target On)
    'targOn',  0.05,   0.2,    0.15;  % 2. Visual (rel. to Target On)
    'pdsFixOff', -0.15,  0.05,   0.2;   % 3. Delay (rel. to Fixation Off)
    'sacOn', -0.025, 0.05,   0.075  % 4. Saccade (rel. to Saccade On)
};

% --- Firing Rate Calculation ---
for i_trial = 1:nTrials
    t_idx = trialInd(i_trial);
    for i_epoch = 1:size(epochs, 1)
        event_name = epochs{i_epoch, 1};
        event_time = eventTimes.(event_name)(t_idx);

        start_offset = epochs{i_epoch, 2};
        end_offset   = epochs{i_epoch, 3};
        win_dur      = epochs{i_epoch, 4};

        if isnan(event_time)
            epoch_frs(:, i_epoch, i_trial) = NaN;
            continue;
        end

        start_time = event_time + start_offset;
        end_time   = event_time + end_offset;

        for i_cluster = 1:nClusters
            cluster_id = cluster_ids(i_cluster);
            spike_times = all_spike_times(all_spike_clusters == cluster_id);
            spike_count = sum(spike_times >= start_time & spike_times < end_time);
            epoch_frs(i_cluster, i_epoch, i_trial) = spike_count / win_dur;
        end
    end
end

% --- Statistical Selection ---
sig_epoch_comparison = false(nClusters, 3); % [Vis vs Base, Delay vs Base, Sac vs Base]
selected_neurons = false(nClusters, 1);

if nTrials < 2
    fprintf('Warning: Too few trials (%d) for statistical selection.\n', nTrials);
    return;
end

for i_cluster = 1:nClusters
    neuron_frs = squeeze(epoch_frs(i_cluster, :, :))';
    neuron_frs = neuron_frs(~any(isnan(neuron_frs), 2), :);

    if size(neuron_frs, 1) < 2
        continue; % Not enough valid trials for this neuron
    end

    try
        [p_friedman, ~] = friedman(neuron_frs, 1, 'off');
        if p_friedman < 0.05
            alpha_corr = 0.05 / 3; % Bonferroni correction
            comparisons = [1 2; 1 3; 1 4];

            for i_comp = 1:size(comparisons, 1)
                col1 = comparisons(i_comp, 1);
                col2 = comparisons(i_comp, 2);
                p_ranksum = ranksum(neuron_frs(:, col1), neuron_frs(:, col2));
                if p_ranksum < alpha_corr
                    sig_epoch_comparison(i_cluster, i_comp) = true;
                end
            end
        end
    catch ME
        fprintf('Stat test failed for cluster %d: %s\n', cluster_ids(i_cluster), ME.message);
    end

    if any(sig_epoch_comparison(i_cluster, :)) && max(mean(neuron_frs, 1, 'omitnan')) > 5
        selected_neurons(i_cluster) = true;
    end
end

% Broaden selection criteria if too few neurons are initially selected
if nnz(selected_neurons) < min(10, nClusters) && nClusters > 0
    fprintf('Too few neurons selected (%d). Broadening criteria...\n', nnz(selected_neurons));
    selected_neurons = any(sig_epoch_comparison, 2);
end

end
