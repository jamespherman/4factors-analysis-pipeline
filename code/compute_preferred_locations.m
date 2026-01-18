%% compute_preferred_locations.m
%
% Computes the preferred target location for each neuron based on visual
% (or saccade) response magnitude. This enables per-neuron trial selection
% for factor comparisons, maximizing sensitivity to factor-related modulation
% by using each neuron's optimal location.
%
% SCIENTIFIC BACKGROUND:
% SC neurons have spatially localized receptive fields (RFs). While the
% experiment aims to place location 1 at the RF center, RF estimation
% isn't always accurate. Some neurons respond best to location 2, 3, or 4.
% Using a pooled contralateral mask dilutes the signal for neurons whose RF
% is at a non-standard location.
%
% For probability: since it's only manipulated at locations 1 & 3 (high
% probability at one, low at the other depending on block), this function
% determines which of these two locations has the stronger visual response.
%
% INPUTS:
%   session_data   - Standard session data structure
%   core_data      - Binned spike data from prepare_core_data
%   conditions     - Trial masks from define_task_conditions
%   condition_defs - Analysis configuration from define_task_conditions
%
% OUTPUTS:
%   preferred_locations - Struct with fields:
%       .visual         - [n_neurons × 1] preferred location based on visual response (1-4)
%       .saccade        - [n_neurons × 1] preferred location based on saccade response (1-4)
%       .for_factors    - [n_neurons × 1] location for reward/salience/identity (1-4)
%       .for_probability- [n_neurons × 1] location for probability (must be 1 or 3)
%       .method_used    - [n_neurons × 1] cell array: 'visual' or 'saccade'
%       .location_responses - [n_neurons × 4] mean response at each location (visual epoch)
%       .location_responses_saccade - [n_neurons × 4] mean response at each location (saccade epoch)
%
% Author: Claude Code
% Date: 2026-01-17

function preferred_locations = compute_preferred_locations(session_data, ...
    core_data, conditions, condition_defs)

%% Setup
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

% Get brain area for window selection
brain_area = session_data.metadata.brain_area;
if ~ismember(brain_area, {'SC', 'SNc'})
    warning('compute_preferred_locations:unknownBrainArea', ...
        'Unknown brain_area ''%s''. Defaulting to SC windows.', brain_area);
    brain_area = 'SC';
end

% Get the epoch windows from window_roc_plan
wroc_plan = condition_defs.window_roc_plan;

% Visual epoch configuration
visual_config = wroc_plan.epoch_windows.visual;
visual_event = visual_config.event;  % 'targetOn'
visual_window = visual_config.(brain_area);

% Saccade epoch configuration
saccade_config = wroc_plan.epoch_windows.perisaccade;
saccade_event = saccade_config.event;  % 'saccadeOnset'
saccade_window = saccade_config.(brain_area);

% Get number of neurons
if ~isfield(core_data.spikes, visual_event)
    error('compute_preferred_locations:missingEvent', ...
        'Visual event ''%s'' not found in core_data.spikes.', visual_event);
end
n_neurons = size(core_data.spikes.(visual_event).rates, 1);

% Get target location indices for each trial
% Note: The conditions struct has been filtered by masterMask in define_task_conditions,
% so we need the targetLocIdx for the same filtered trials.
% We can derive this by finding the 4factors trial indices
codes = initCodes();
isGSac4factors = session_data.trialInfo.taskCode == codes.uniqueTaskCode_gSac_4factors;
isSuccessful = ~cellfun(@isempty, session_data.eventTimes.rewardCell);
masterMask = isGSac4factors & isSuccessful;

% Get targetLocIdx filtered by masterMask (same as conditions)
targetLocIdx = session_data.trialInfo.targetLocIdx(masterMask);

%% Compute Location Responses - Visual Epoch
visual_rates = core_data.spikes.(visual_event).rates;  % [n_neurons × n_trials × n_time_bins]
visual_time = core_data.spikes.(visual_event).time_vector;

% Find bins within the visual analysis window
visual_bin_mask = (visual_time >= visual_window(1)) & (visual_time <= visual_window(2));

if sum(visual_bin_mask) == 0
    error('compute_preferred_locations:noVisualBins', ...
        'No bins found within visual window [%.3f, %.3f].', ...
        visual_window(1), visual_window(2));
end

% Average firing rate in visual window: [n_neurons × n_trials]
visual_window_rates = mean(visual_rates(:, :, visual_bin_mask), 3, 'omitnan');

% Compute mean response at each location (1-4)
location_responses_visual = NaN(n_neurons, 4);
for loc = 1:4
    loc_mask = (targetLocIdx == loc);
    if sum(loc_mask) > 0
        location_responses_visual(:, loc) = mean(visual_window_rates(:, loc_mask), 2, 'omitnan');
    end
end

%% Compute Location Responses - Saccade Epoch
if isfield(core_data.spikes, saccade_event)
    saccade_rates = core_data.spikes.(saccade_event).rates;
    saccade_time = core_data.spikes.(saccade_event).time_vector;

    % Find bins within the saccade analysis window
    saccade_bin_mask = (saccade_time >= saccade_window(1)) & (saccade_time <= saccade_window(2));

    if sum(saccade_bin_mask) > 0
        % Average firing rate in saccade window: [n_neurons × n_trials]
        saccade_window_rates = mean(saccade_rates(:, :, saccade_bin_mask), 3, 'omitnan');

        % Compute mean response at each location (1-4)
        location_responses_saccade = NaN(n_neurons, 4);
        for loc = 1:4
            loc_mask = (targetLocIdx == loc);
            if sum(loc_mask) > 0
                location_responses_saccade(:, loc) = mean(saccade_window_rates(:, loc_mask), 2, 'omitnan');
            end
        end
    else
        location_responses_saccade = NaN(n_neurons, 4);
    end
else
    location_responses_saccade = NaN(n_neurons, 4);
end

%% Determine Preferred Locations
% Threshold for "don't differentiate": max < 1.2 × min (less than 20% difference)
DIFFERENTIATION_RATIO = 1.2;

% Initialize outputs
preferred_visual = NaN(n_neurons, 1);
preferred_saccade = NaN(n_neurons, 1);
preferred_for_factors = NaN(n_neurons, 1);
preferred_for_probability = NaN(n_neurons, 1);
method_used = cell(n_neurons, 1);

for i_neuron = 1:n_neurons
    % Get visual responses for this neuron
    vis_resp = location_responses_visual(i_neuron, :);
    sac_resp = location_responses_saccade(i_neuron, :);

    % Find location with max visual response
    [max_vis, preferred_visual(i_neuron)] = max(vis_resp);

    % Find location with max saccade response
    [~, preferred_saccade(i_neuron)] = max(sac_resp);

    % Check if visual responses differentiate locations
    valid_vis = vis_resp(~isnan(vis_resp));
    if ~isempty(valid_vis)
        min_vis = min(valid_vis);
        % Avoid division by zero: if min is 0 but max > 0, they differentiate
        if min_vis > 0
            vis_differentiates = (max_vis / min_vis) >= DIFFERENTIATION_RATIO;
        else
            vis_differentiates = max_vis > 0;  % If min=0 and max>0, they differentiate
        end
    else
        vis_differentiates = false;
    end

    % Determine method and preferred location for factors
    if vis_differentiates
        preferred_for_factors(i_neuron) = preferred_visual(i_neuron);
        method_used{i_neuron} = 'visual';
    else
        % Fall back to saccade response
        valid_sac = sac_resp(~isnan(sac_resp));
        if ~isempty(valid_sac)
            [max_sac, ~] = max(valid_sac);
            min_sac = min(valid_sac);
            if min_sac > 0
                sac_differentiates = (max_sac / min_sac) >= DIFFERENTIATION_RATIO;
            else
                sac_differentiates = max_sac > 0;
            end

            if sac_differentiates
                preferred_for_factors(i_neuron) = preferred_saccade(i_neuron);
                method_used{i_neuron} = 'saccade';
            else
                % Neither differentiates - default to visual max anyway
                preferred_for_factors(i_neuron) = preferred_visual(i_neuron);
                method_used{i_neuron} = 'visual';
            end
        else
            % No valid saccade data - use visual
            preferred_for_factors(i_neuron) = preferred_visual(i_neuron);
            method_used{i_neuron} = 'visual';
        end
    end

    % Determine probability location (must be 1 or 3)
    pref_loc = preferred_for_factors(i_neuron);
    if pref_loc == 1 || pref_loc == 3
        % Use the preferred location directly
        preferred_for_probability(i_neuron) = pref_loc;
    else
        % Preferred location is 2 or 4 - compare responses at 1 vs 3
        resp_at_1 = vis_resp(1);
        resp_at_3 = vis_resp(3);

        if isnan(resp_at_1) && isnan(resp_at_3)
            % No data at either location - default to 1
            preferred_for_probability(i_neuron) = 1;
        elseif isnan(resp_at_1)
            preferred_for_probability(i_neuron) = 3;
        elseif isnan(resp_at_3)
            preferred_for_probability(i_neuron) = 1;
        elseif resp_at_1 >= resp_at_3
            preferred_for_probability(i_neuron) = 1;
        else
            preferred_for_probability(i_neuron) = 3;
        end
    end

    % Handle edge case: NaN preferred location (shouldn't happen, but be safe)
    if isnan(preferred_for_factors(i_neuron))
        preferred_for_factors(i_neuron) = 1;  % Default to location 1
        method_used{i_neuron} = 'visual';
    end
    if isnan(preferred_for_probability(i_neuron))
        preferred_for_probability(i_neuron) = 1;
    end
end

%% Package Output
preferred_locations = struct();
preferred_locations.visual = preferred_visual;
preferred_locations.saccade = preferred_saccade;
preferred_locations.for_factors = preferred_for_factors;
preferred_locations.for_probability = preferred_for_probability;
preferred_locations.method_used = method_used;
preferred_locations.location_responses = location_responses_visual;
preferred_locations.location_responses_saccade = location_responses_saccade;

% Add metadata about the computation
preferred_locations.metadata = struct();
preferred_locations.metadata.brain_area = brain_area;
preferred_locations.metadata.visual_window = visual_window;
preferred_locations.metadata.saccade_window = saccade_window;
preferred_locations.metadata.differentiation_ratio = DIFFERENTIATION_RATIO;

end
