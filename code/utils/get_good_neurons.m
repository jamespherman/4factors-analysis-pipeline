function good_neuron_ids = get_good_neurons(session_data)
% GET_GOOD_NEURONS Filters for neurons marked as 'Good' quality.
%
%   Usage:
%   good_neuron_ids = get_good_neurons(session_data)
%
%   Inputs:
%   session_data - A struct containing session data, which must include
%                  session_data.spikes.cluster_info.
%
%   Outputs:
%   good_neuron_ids - A numeric column vector of 'Good' neuron IDs.

if isfield(session_data, 'spikes') && isfield(session_data.spikes, 'cluster_info')
    cluster_info = session_data.spikes.cluster_info;

    % Find rows where the 'group' is 'Good'
    good_neurons_mask = strcmp(cluster_info.group, 'Good');

    % Get the corresponding cluster IDs
    good_neuron_ids = cluster_info.cluster_id(good_neurons_mask);
else
    % Return empty if the required data is not present
    good_neuron_ids = [];
    warning('Could not find cluster_info to identify good neurons.');
end

% Ensure the output is a column vector
if ~iscolumn(good_neuron_ids)
    good_neuron_ids = good_neuron_ids';
end

end
