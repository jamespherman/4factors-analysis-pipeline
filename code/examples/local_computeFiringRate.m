function firingRate = local_computeFiringRate(sc, startTime, endTime, duration, nNeurons)
    % computes firing rate for all neurons within a specified time window.
    if nNeurons == 0
        firingRate = zeros(0,1); % return correctly dimensioned empty array
        return;
    end
    firingRate = zeros(nNeurons, 1); % initialize

    if isempty(sc.clusterTimes) % no spike times available
        return;
    end

    % iterate through each neuron
    for i_neuron = 1:nNeurons
        % check if current neuron index is valid and has spike times
        if i_neuron <= numel(sc.clusterTimes) && ~isempty(sc.clusterTimes{i_neuron})
            neuron_spikes_in_window = sum(sc.clusterTimes{i_neuron} >= startTime & sc.clusterTimes{i_neuron} <= endTime);
            if duration > 0
                firingRate(i_neuron) = neuron_spikes_in_window / duration;
            else
                firingRate(i_neuron) = 0; % avoid division by zero
            end
        else
            firingRate(i_neuron) = 0; % no spike times for this neuron or index out of bounds
        end
    end
end