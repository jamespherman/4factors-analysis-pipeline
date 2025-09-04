function selected_neuron_ids = select_neurons(session_data)
% SELECT_NEURONS acts as a controller to select neurons for analysis.
%
%   This function is a wrapper that calls area-specific screening
%   functions based on the brain_area specified in the session's
%   metadata.
%
%   Usage:
%   selected_neuron_ids = select_neurons(session_data)
%
%   Inputs:
%   session_data - A struct containing session data, including metadata.
%                  Must contain session_data.metadata.brain_area.
%
%   Outputs:
%   selected_neuron_ids - A numeric column vector of selected neuron IDs.

% Get the full list of neuron IDs from the session data.
all_neuron_ids = session_data.spikes.cluster_info.cluster_id;

% Initialize the mask to select no neurons by default.
n_neurons = numel(all_neuron_ids);
is_selected_mask = false(n_neurons, 1);

% Use a switch statement to call the appropriate screening function
% based on the brain area specified in the metadata.
switch session_data.metadata.brain_area
    case 'SNc'
        % For SNc, screen for putative dopamine neurons.
        is_selected_mask = analysis.screen_da_neurons(session_data, ...
            session_data.metadata.unique_id);

    case 'SC'
        % For SC, screen for task-modulated neurons.
        % The hemisphere is needed to determine contralateral targets.
        if isfield(session_data.metadata, 'hemisphere')
            hemisphere = session_data.metadata.hemisphere;
            is_selected_mask = analysis.screen_sc_neurons(session_data, ...
                hemisphere);
        elseif isfield(session_data.metadata, 'scSide') % Fallback check
            hemisphere = session_data.metadata.scSide;
            is_selected_mask = analysis.screen_sc_neurons(session_data, ...
                hemisphere);
        else
            % If hemisphere is not specified, we cannot select neurons.
            warning(['No hemisphere or scSide field found in metadata ' ...
                     'for SC area. Returning empty list.']);
        end

    otherwise
        % If the brain area is not recognized, return an empty list.
        warning(['No specific selection criteria for brain area: %s. ' ...
                 'Returning empty list.'], ...
                 session_data.metadata.brain_area);
end

% Use the logical mask to get the final list of neuron IDs.
selected_neuron_ids = all_neuron_ids(is_selected_mask);

% Ensure the output is a column vector, just in case.
if ~iscolumn(selected_neuron_ids)
    selected_neuron_ids = selected_neuron_ids';
end

end
