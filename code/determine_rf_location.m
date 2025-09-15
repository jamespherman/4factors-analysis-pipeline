function [rf_location_idx, hemisphere] = determineRfLocation(session_data)
% determineRfLocation Programmatically determines the target location
% corresponding to the SC receptive field (RF) for a given session.
%
%   Note: The function name `determineRfLocation` uses camelCase to adhere
%   to the project's coding standards (see AGENTS.md), rather than the
%   snake_case originally requested.
%
%   Usage:
%   [rf_location_idx, hemisphere] = determineRfLocation(session_data)
%
%   Inputs:
%   session_data - A struct containing session-specific data. If the
%                  brain_area is 'SNc', the function will automatically
%                  find and use the corresponding 'SC' session data.
%
%   Outputs:
%   rf_location_idx - The numeric index (1-4) of the gSac_4factors target
%                     location that falls within the receptive field.
%   hemisphere      - A string ('left' or 'right') indicating the
%                     hemisphere of the recorded SC.

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', ...
    '4factors-analysis-pipeline');
addpath(fullfile(project_root, 'code', 'utils'));

%% Initialize outputs
rf_location_idx = NaN;
hemisphere = 'unknown';

%% Handle SNc sessions by finding the SC partner
% If the provided session is an SNc recording, we need to find its
% simultaneously recorded SC partner to determine the RF, as the RF is
% defined by the SC activity.

% Make a copy of the input session_data to avoid modifying the original
sc_session_data = session_data;

if strcmp(sc_session_data.metadata.brain_area, 'SNc')
    fprintf(['Input session is SNc. Finding corresponding SC ' ...
        'session...\n']);

    % Define paths for manifest and data
    manifest_path = fullfile(project_root, 'config', ...
        'session_manifest.csv');
    data_base_path = fullfile(project_root, 'data');

    if ~isfile(manifest_path)
        error('determineRfLocation:ManifestNotFound', ...
            'Could not find session_manifest.csv at: %s', manifest_path);
    end

    % Read the manifest
    manifest = readtable(manifest_path);

    % Find the session_group_id for the current SNc session
    current_id = sc_session_data.metadata.unique_id;
    senc_row = manifest(strcmp(manifest.unique_id, current_id), :);
    if isempty(senc_row)
        error('determineRfLocation:IdNotFound', ...
            'Current unique_id "%s" not found in manifest.', current_id);
    end
    session_group_id = senc_row.session_group_id{1};

    % Find the partner SC session in the same group
    sc_partner_row = manifest(...
        strcmp(manifest.session_group_id, session_group_id) & ...
        strcmp(manifest.brain_area, 'SC'), :);

    if isempty(sc_partner_row)
        warning(['No corresponding SC session found for group %s. ' ...
            'Cannot determine RF location.'], session_group_id);
        return; % Return with default NaN/unknown values
    elseif height(sc_partner_row) > 1
        warning(['Multiple SC sessions found for group %s. Using the ' ...
            'first one.'], session_group_id);
        sc_partner_row = sc_partner_row(1, :);
    end

    sc_unique_id = sc_partner_row.unique_id{1};
    fprintf('Found SC partner session: %s\n', sc_unique_id);

    % Load the SC partner session data
    sc_session_data = data_handling.load_session(sc_unique_id, ...
        manifest_path, data_base_path);
end

%% Core Logic: Identify RF from SC activity
% Call screen_sc_neurons to get statistical results for each trial group.
[~, ~, hemisphere, trial_groups, group_sig_results] = ...
    analysis.screen_sc_neurons(sc_session_data);

% Find trial groups corresponding to the gSac_4factors task.
codes = initCodes();
is_4factors_trial = sc_session_data.trialInfo.taskCode == ...
    codes.uniqueTaskCode_gSac_4factors;

% Get the target locations (thetas) for all 4factors trials
thetas_4factors = sc_session_data.trialInfo.targetTheta(is_4factors_trial);
unique_thetas = unique(thetas_4factors);
n_locations = length(unique_thetas);

if n_locations == 0
    warning('No gSac_4factors trials found in this session.');
    return;
end

% Thetas are not guaranteed to be sorted, so sort them to ensure a
% consistent index mapping (1 -> smallest theta, 4 -> largest theta).
sorted_thetas = sort(unique_thetas);

% Calculate a "modulation score" for each target location.
modulation_scores = zeros(1, n_locations);
thetas_all = sc_session_data.trialInfo.targetTheta;

for i_group = 1:length(trial_groups)
    group_mask = trial_groups{i_group};

    % Check if this group belongs to the 4factors task by seeing if the
    % first trial in the group is a 4factors trial.
    first_trial_idx = find(group_mask, 1);
    if isempty(first_trial_idx) || ~is_4factors_trial(first_trial_idx)
        continue; % Skip non-4factors groups (e.g., gSac_jph)
    end

    % Determine the theta for this group
    group_theta = unique(thetas_all(group_mask));
    if isempty(group_theta) || isnan(group_theta)
        continue;
    end

    % Find the index of this theta in our sorted list
    [~, loc_idx] = ismember(group_theta, sorted_thetas);
    if loc_idx == 0
        continue; % Should not happen, but as a safeguard
    end

    % The score is the total number of significant epoch modulations for
    % this location across all neurons.
    modulation_scores(loc_idx) = sum(group_sig_results{i_group}, 'all');
end

% The RF location is the one with the highest modulation score.
[~, winning_idx] = max(modulation_scores);

if isempty(winning_idx)
    warning('Could not determine a winning RF location.');
    return;
end

rf_location_idx = winning_idx;
fprintf('Determined RF location index: %d\n', rf_location_idx);

end
