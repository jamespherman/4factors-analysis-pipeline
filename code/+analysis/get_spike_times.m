function spike_times_cell = get_spike_times(session_data, neuron_ids)
% GET_SPIKE_TIMES Extracts spike times for a specified list of neurons.
%
%   Usage:
%   spike_times_cell = get_spike_times(session_data, neuron_ids)
%
%   Inputs:
%   session_data - The main data struct for a session. It must contain
%                  a 'spikes' field with 'clusters' and 'times'.
%   neuron_ids   - A numeric vector of the cluster_ids for which to
%                  extract spike times.
%
%   Outputs:
%   spike_times_cell - A cell array where each cell contains the spike
%                      time vector for one neuron. The order of cells
%                      matches the order of the input neuron_ids.

% Initialize an empty cell array with the same size as neuron_ids.
spike_times_cell = cell(size(neuron_ids));

% Loop through each id in the input neuron_ids vector.
for i = 1:length(neuron_ids)
    current_id = neuron_ids(i);
    
    % Create a logical mask to find all spikes belonging to that cluster.
    mask = session_data.spikes.clusters == current_id;
    
    % Use this mask to select the corresponding timestamps.
    spike_times = session_data.spikes.times(mask);
    
    % Store the resulting vector of spike times in the cell array.
    spike_times_cell{i} = spike_times;
end

end
