function scSide = local_defineScSide(sc, gSac, giveFeed)
    giveFeed('local_defineScSide: Determining SC recording side...');
    % define window for assessing visual response to target onset
    duration = 0.15; % visual window duration: targetOn + 0.05 to targetOn + 0.20
    nNeurons = sc.nClusters;

    if nNeurons == 0
        giveFeed('WARNING in local_defineScSide: No clusters found. Defaulting scSide to "unknown".');
        scSide = 'unknown';
        return;
    end

    % identify high saliency trials for left (T1) and right (T2) targets
    highSalTrials = gSac.isHighSaliencyTrial;
    T1_idx = find(highSalTrials & gSac.isT1Trial); % typically Left Target
    T2_idx = find(highSalTrials & gSac.isT2Trial); % typically Right Target

    if isempty(T1_idx) || isempty(T2_idx)
        giveFeed('WARNING in local_defineScSide: Not enough T1 or T2 high-saliency trials to determine SC side. Defaulting to "right".');
        scSide = 'right';
        return;
    end

    fr_T1 = zeros(numel(T1_idx), 1); % to store mean firing rates for T1 trials
    fr_T2 = zeros(numel(T2_idx), 1); % to store mean firing rates for T2 trials

    % calculate mean firing rate across neurons for T1 trials
    for i_trial = 1:numel(T1_idx)
        t = T1_idx(i_trial);
        startTime = gSac.targetOnTime(t) + 0.05;
        endTime = gSac.targetOnTime(t) + 0.20;
        
        trial_spike_counts_per_neuron = zeros(nNeurons, 1);
        if ~isempty(sc.clusterTimes) && nNeurons > 0
            for i_neuron = 1:nNeurons
                if i_neuron <= numel(sc.clusterTimes) && ~isempty(sc.clusterTimes{i_neuron})
                     neuron_spikes_in_window = sc.clusterTimes{i_neuron}(sc.clusterTimes{i_neuron} >= startTime & sc.clusterTimes{i_neuron} <= endTime);
                     trial_spike_counts_per_neuron(i_neuron) = numel(neuron_spikes_in_window);
                end
            end
        end
        fr_T1(i_trial) = mean(trial_spike_counts_per_neuron) / duration; % mean FR across all neurons for this T1 trial
    end

    % calculate mean firing rate across neurons for T2 trials
    for i_trial = 1:numel(T2_idx)
        t = T2_idx(i_trial);
        startTime = gSac.targetOnTime(t) + 0.05;
        endTime = gSac.targetOnTime(t) + 0.20;
        
        trial_spike_counts_per_neuron = zeros(nNeurons, 1);
         if ~isempty(sc.clusterTimes) && nNeurons > 0
            for i_neuron = 1:nNeurons
                if i_neuron <= numel(sc.clusterTimes) && ~isempty(sc.clusterTimes{i_neuron})
                    neuron_spikes_in_window = sc.clusterTimes{i_neuron}(sc.clusterTimes{i_neuron} >= startTime & sc.clusterTimes{i_neuron} <= endTime);
                    trial_spike_counts_per_neuron(i_neuron) = numel(neuron_spikes_in_window);
                end
            end
        end
        fr_T2(i_trial) = mean(trial_spike_counts_per_neuron) / duration; % mean FR across all neurons for this T2 trial
    end

    % determine SC side based on higher mean firing rate
    mean_fr_T1 = mean(fr_T1(~isnan(fr_T1)));
    mean_fr_T2 = mean(fr_T2(~isnan(fr_T2)));

    if isnan(mean_fr_T1) && isnan(mean_fr_T2)
        giveFeed('WARNING in local_defineScSide: Mean FR for T1 and T2 are NaN. Defaulting to "right".');
        scSide = 'right';
        return;
    elseif isnan(mean_fr_T1)
        scSide = 'left'; % if T1 FR is NaN, assume T2 response was higher (Right SC -> Left Target/T1; Left SC -> Right Target/T2)
        giveFeed('WARNING in local_defineScSide: Mean FR T1 is NaN. Defaulting SC side to left (implies T2 response dominant).');
        return;
    elseif isnan(mean_fr_T2)
        scSide = 'right'; % if T2 FR is NaN, assume T1 response was higher
        giveFeed('WARNING in local_defineScSide: Mean FR T2 is NaN. Defaulting SC side to right (implies T1 response dominant).');
        return;
    end

    if mean_fr_T1 > mean_fr_T2 % T1 (Left Target) response is stronger
        scSide = 'right'; % right SC is contralateral to left targets
    else % T2 (Right Target) response is stronger or equal
        scSide = 'left';  % left SC is contralateral to right targets
    end
    giveFeed(['local_defineScSide: Mean FR T1 (Left Target): ' num2str(mean_fr_T1) ', Mean FR T2 (Right Target): ' num2str(mean_fr_T2)]);
end