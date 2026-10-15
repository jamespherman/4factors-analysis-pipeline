function plot_sc_spatial_selectivity()
% plot_sc_spatial_selectivity - Plots spatial selectivity of Classic vs Non-classic SC neurons.
%
% Generates two figures:
%   1. Contra/Ipsi Collapsed (116x116 px): Grouped by neuron type, then
%      contra/ipsi, with activated/suppressed bars.
%   2. All Four Locations (200x116 px): Grouped by neuron type, then
%      location (1-4), with activated/suppressed bars.
%
% Uses Wilson score 95% CIs for error bars.
%
% Output files:
%   figures/sc_spatial_selectivity_collapsed.pdf
%   figures/sc_spatial_selectivity_4locations.pdf

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
addpath(fullfile(project_root, 'code', 'utils'));

%% Load Aggregated Data for Session List
fprintf('Loading aggregated analysis data...\n');
agg_data_path = fullfile(project_root, 'data', 'processed', 'aggregated_analysis_data.mat');
if ~exist(agg_data_path, 'file')
    error('Aggregated data not found. Run aggregate_analysis_results first.');
end
load(agg_data_path, 'aggregated_sc_data', 'session_ids');

% Get neuron metrics table for neuron types
nmt = aggregated_sc_data.neuron_metrics_table;

%% Define Standard Angles
% Location 1 = 45°, Location 2 = 135°, Location 3 = 225°, Location 4 = 315°
standard_angles = [45, 135, 225, 315];

%% Initialize Counters
% We need to track per RF-relative location:
%   - Location 1 = RF center (upper contralateral)
%   - Location 2 = lower contralateral
%   - Location 3 = upper ipsilateral
%   - Location 4 = lower ipsilateral
%
% For collapsed version:
%   - Contra = Loc 1 + Loc 2
%   - Ipsi = Loc 3 + Loc 4

% Counts structure: counts.(neuron_type).(location_type).(response_type)
% Initialize all counters
neuron_types = {'classic', 'interneuron'};
rf_locations = {'loc1', 'loc2', 'loc3', 'loc4'};
contra_ipsi = {'contra', 'ipsi'};
response_types = {'activated', 'suppressed', 'total'};

% For 4-location analysis
loc_counts = struct();
for nt = 1:length(neuron_types)
    for loc = 1:4
        loc_name = sprintf('loc%d', loc);
        loc_counts.(neuron_types{nt}).(loc_name).activated = 0;
        loc_counts.(neuron_types{nt}).(loc_name).suppressed = 0;
        loc_counts.(neuron_types{nt}).(loc_name).total = 0;
    end
end

% For collapsed contra/ipsi analysis
ci_counts = struct();
for nt = 1:length(neuron_types)
    for ci = 1:length(contra_ipsi)
        ci_counts.(neuron_types{nt}).(contra_ipsi{ci}).activated = 0;
        ci_counts.(neuron_types{nt}).(contra_ipsi{ci}).suppressed = 0;
        ci_counts.(neuron_types{nt}).(contra_ipsi{ci}).total = 0;
    end
end

%% Process Each Session
fprintf('Processing sessions...\n');
for i_sess = 1:length(session_ids.sc)
    session_id = session_ids.sc{i_sess};

    % Load session data
    session_data_path = fullfile(project_root, 'data', 'processed', ...
        session_id, [session_id '_session_data.mat']);

    if ~exist(session_data_path, 'file')
        fprintf('  Skipping %s (file not found)\n', session_id);
        continue;
    end

    data = load(session_data_path, 'session_data');
    sd = data.session_data;

    % Check for required fields
    if ~isfield(sd, 'cached_friedman_results') || isempty(sd.cached_friedman_results)
        fprintf('  Skipping %s (no cached_friedman_results)\n', session_id);
        continue;
    end

    % Get scSide
    if isfield(sd, 'analysis') && isfield(sd.analysis, 'scSide')
        scSide = sd.analysis.scSide;
    else
        fprintf('  Skipping %s (no scSide)\n', session_id);
        continue;
    end

    friedman = sd.cached_friedman_results;

    % Get neuron classifications for this session
    session_mask = strcmp(nmt.session_id, session_id);
    session_neuron_types = nmt.neuron_type(session_mask);

    % Validate dimensions
    n_clusters = size(friedman.sig_by_location, 1);
    if length(session_neuron_types) ~= n_clusters
        fprintf('  Warning: dimension mismatch for %s\n', session_id);
        continue;
    end

    % Define contralateral angles based on scSide
    % Right SC: contralateral = left visual field = 135° and 225°
    % Left SC: contralateral = right visual field = 45° and 315°
    if strcmp(scSide, 'right')
        contra_angles = [135, 225];
        ipsi_angles = [45, 315];
    elseif strcmp(scSide, 'left')
        contra_angles = [45, 315];
        ipsi_angles = [135, 225];
    else
        fprintf('  Skipping %s (unknown scSide: %s)\n', session_id, scSide);
        continue;
    end

    % Map absolute angles to RF-relative locations
    % RF-relative mapping:
    %   Loc 1 = upper contralateral (RF center)
    %   Loc 2 = lower contralateral
    %   Loc 3 = upper ipsilateral
    %   Loc 4 = lower ipsilateral
    %
    % For right SC (contra = left visual field):
    %   Loc 1 = 135° (upper left), Loc 2 = 225° (lower left)
    %   Loc 3 = 45° (upper right), Loc 4 = 315° (lower right)
    %
    % For left SC (contra = right visual field):
    %   Loc 1 = 45° (upper right), Loc 2 = 315° (lower right)
    %   Loc 3 = 135° (upper left), Loc 4 = 225° (lower left)

    if strcmp(scSide, 'right')
        angle_to_rf_loc = containers.Map([135, 225, 45, 315], [1, 2, 3, 4]);
    else % left
        angle_to_rf_loc = containers.Map([45, 315, 135, 225], [1, 2, 3, 4]);
    end

    % Location angles from friedman results
    location_angles = friedman.location_angles;
    n_groups = length(location_angles);

    fprintf('  Processing %s (scSide=%s, %d neurons)\n', session_id, scSide, n_clusters);

    % Process each neuron
    for i_neuron = 1:n_clusters
        neuron_type = session_neuron_types{i_neuron};

        % Skip excluded neurons
        if startsWith(neuron_type, 'excluded')
            continue;
        end

        % Map neuron type
        if strcmp(neuron_type, 'classic')
            nt_key = 'classic';
        elseif strcmp(neuron_type, 'interneuron')
            nt_key = 'interneuron';
        else
            continue;
        end

        % Process each location group
        for i_group = 1:n_groups
            angle = location_angles(i_group);

            % Map angle to standard location index
            [~, std_loc_idx] = min(abs(standard_angles - angle));

            % Skip non-standard angles (e.g., jph trials)
            if abs(standard_angles(std_loc_idx) - angle) > 5
                continue;
            end

            % Get RF-relative location
            actual_angle = standard_angles(std_loc_idx);
            if ~isKey(angle_to_rf_loc, actual_angle)
                continue;
            end
            rf_loc = angle_to_rf_loc(actual_angle);
            rf_loc_name = sprintf('loc%d', rf_loc);

            % Determine contra/ipsi
            if ismember(actual_angle, contra_angles)
                ci_key = 'contra';
            else
                ci_key = 'ipsi';
            end

            % Check for activation or suppression (visual epoch = index 1)
            has_activation = friedman.sig_by_location(i_neuron, i_group, 1);
            has_suppression = false;
            if isfield(friedman, 'suppression_by_location')
                has_suppression = friedman.suppression_by_location(i_neuron, i_group, 1);
            end

            % Update location counts
            loc_counts.(nt_key).(rf_loc_name).total = loc_counts.(nt_key).(rf_loc_name).total + 1;
            if has_activation
                loc_counts.(nt_key).(rf_loc_name).activated = loc_counts.(nt_key).(rf_loc_name).activated + 1;
            end
            if has_suppression
                loc_counts.(nt_key).(rf_loc_name).suppressed = loc_counts.(nt_key).(rf_loc_name).suppressed + 1;
            end

            % Update collapsed counts
            ci_counts.(nt_key).(ci_key).total = ci_counts.(nt_key).(ci_key).total + 1;
            if has_activation
                ci_counts.(nt_key).(ci_key).activated = ci_counts.(nt_key).(ci_key).activated + 1;
            end
            if has_suppression
                ci_counts.(nt_key).(ci_key).suppressed = ci_counts.(nt_key).(ci_key).suppressed + 1;
            end
        end
    end
end

%% Print Summary
fprintf('\n=== Spatial Selectivity Summary ===\n\n');

fprintf('4-Location Counts (RF-relative):\n');
for nt = 1:length(neuron_types)
    fprintf('  %s:\n', neuron_types{nt});
    for loc = 1:4
        loc_name = sprintf('loc%d', loc);
        lc = loc_counts.(neuron_types{nt}).(loc_name);
        fprintf('    Loc %d: %d/%d activated (%.1f%%), %d/%d suppressed (%.1f%%)\n', ...
            loc, lc.activated, lc.total, 100*lc.activated/max(lc.total,1), ...
            lc.suppressed, lc.total, 100*lc.suppressed/max(lc.total,1));
    end
end

fprintf('\nCollapsed Counts:\n');
for nt = 1:length(neuron_types)
    fprintf('  %s:\n', neuron_types{nt});
    for ci = 1:length(contra_ipsi)
        cc = ci_counts.(neuron_types{nt}).(contra_ipsi{ci});
        fprintf('    %s: %d/%d activated (%.1f%%), %d/%d suppressed (%.1f%%)\n', ...
            contra_ipsi{ci}, cc.activated, cc.total, 100*cc.activated/max(cc.total,1), ...
            cc.suppressed, cc.total, 100*cc.suppressed/max(cc.total,1));
    end
end

%% Prepare Data for Plotting

% Wilson score CI function
wilson_ci = @(k, n) wilson_score_ci(k, n, 0.05);

% Color scheme
color_activated = [0.3 0.3 0.3];  % Dark gray for activated
color_suppressed = [0.7 0.7 0.7]; % Light gray for suppressed

%% Figure 1: Contra/Ipsi Collapsed (116 x 116 pixels)
fprintf('\nGenerating collapsed figure...\n');

% Bar order: Classic-Contra-Act, Classic-Contra-Sup, Classic-Ipsi-Act, Classic-Ipsi-Sup,
%            NonClassic-Contra-Act, NonClassic-Contra-Sup, NonClassic-Ipsi-Act, NonClassic-Ipsi-Sup

% Collect data
x_positions = [];
y_values = [];
y_ci = [];
bar_colors = [];

bar_idx = 0;
group_spacing = 1.0;  % Space between Classic and Non-classic groups
subgroup_spacing = 0.5;  % Space between Contra and Ipsi
bar_spacing = 0.3;  % Space between Act and Sup bars

for nt = 1:length(neuron_types)
    for ci = 1:length(contra_ipsi)
        for resp = 1:2  % 1=activated, 2=suppressed
            bar_idx = bar_idx + 1;

            cc = ci_counts.(neuron_types{nt}).(contra_ipsi{ci});

            if resp == 1
                k = cc.activated;
                bar_color = color_activated;
            else
                k = cc.suppressed;
                bar_color = color_suppressed;
            end
            n = cc.total;

            % Calculate proportion and CI
            if n > 0
                p = k / n;
                [ci_lo, ci_hi] = wilson_ci(k, n);
            else
                p = 0;
                ci_lo = 0;
                ci_hi = 0;
            end

            % Calculate x position
            x = (nt - 1) * (2 * subgroup_spacing + 2 * bar_spacing + group_spacing) + ...
                (ci - 1) * (bar_spacing + subgroup_spacing) + ...
                (resp - 1) * bar_spacing;

            x_positions(bar_idx) = x;
            y_values(bar_idx) = p;
            y_ci(bar_idx, :) = [ci_lo, ci_hi];
            bar_colors(bar_idx, :) = bar_color;
        end
    end
end

% Create figure
fig1 = figure('Color', 'w', 'Position', [100, 100, 116, 116]);
ax1 = axes('Parent', fig1);
hold(ax1, 'on');

% Plot bars using mybarerr
for i = 1:length(x_positions)
    mybarerr(x_positions(i), y_values(i), y_ci(i,:), [], bar_colors(i,:), 0.12);
end

% Set axis limits
x_min = min(x_positions) - 0.3;
x_max = max(x_positions) + 0.3;
y_min = 0;
y_max = 1;

xlim(ax1, [x_min, x_max]);
ylim(ax1, [y_min, y_max]);

% Draw disconnected axes
X_axis = [x_min - 0.05*(x_max-x_min), x_min, x_max];
Y_axis = [y_min - 0.05*(y_max-y_min), y_min, y_max];
drawAxes(ax1, X_axis, Y_axis, 0.04, 0, 3);

hold(ax1, 'off');

% Save figure
output_dir = fullfile(project_root, 'figures');
fig1_path = fullfile(output_dir, 'sc_spatial_selectivity_collapsed.pdf');
pdfSave(fig1_path, fig1.Position(3:4)/72, fig1);
fprintf('Saved: %s\n', fig1_path);
fprintf('Y-axis limits: [%.2f, %.2f]\n', y_min, y_max);

close(fig1);

%% Figure 2: All Four Locations (200 x 116 pixels)
fprintf('\nGenerating 4-location figure...\n');

% Bar order: For each neuron type (Classic, Non-classic):
%   For each location (1, 2, 3, 4):
%     Activated, Suppressed

x_positions_4loc = [];
y_values_4loc = [];
y_ci_4loc = [];
bar_colors_4loc = [];

bar_idx = 0;
group_spacing = 1.2;    % Space between Classic and Non-classic
loc_spacing = 0.5;      % Space between locations
bar_spacing = 0.25;     % Space between Act and Sup bars

for nt = 1:length(neuron_types)
    for loc = 1:4
        loc_name = sprintf('loc%d', loc);
        for resp = 1:2
            bar_idx = bar_idx + 1;

            lc = loc_counts.(neuron_types{nt}).(loc_name);

            if resp == 1
                k = lc.activated;
                bar_color = color_activated;
            else
                k = lc.suppressed;
                bar_color = color_suppressed;
            end
            n = lc.total;

            if n > 0
                p = k / n;
                [ci_lo, ci_hi] = wilson_ci(k, n);
            else
                p = 0;
                ci_lo = 0;
                ci_hi = 0;
            end

            % Calculate x position
            x = (nt - 1) * (4 * loc_spacing + 4 * bar_spacing + group_spacing) + ...
                (loc - 1) * (bar_spacing + loc_spacing) + ...
                (resp - 1) * bar_spacing;

            x_positions_4loc(bar_idx) = x;
            y_values_4loc(bar_idx) = p;
            y_ci_4loc(bar_idx, :) = [ci_lo, ci_hi];
            bar_colors_4loc(bar_idx, :) = bar_color;
        end
    end
end

% Create figure
fig2 = figure('Color', 'w', 'Position', [100, 100, 200, 116]);
ax2 = axes('Parent', fig2);
hold(ax2, 'on');

% Plot bars
for i = 1:length(x_positions_4loc)
    mybarerr(x_positions_4loc(i), y_values_4loc(i), y_ci_4loc(i,:), [], bar_colors_4loc(i,:), 0.1);
end

% Set axis limits
x_min_4 = min(x_positions_4loc) - 0.3;
x_max_4 = max(x_positions_4loc) + 0.3;
y_min_4 = 0;
y_max_4 = 1;

xlim(ax2, [x_min_4, x_max_4]);
ylim(ax2, [y_min_4, y_max_4]);

% Draw disconnected axes
X_axis_4 = [x_min_4 - 0.05*(x_max_4-x_min_4), x_min_4, x_max_4];
Y_axis_4 = [y_min_4 - 0.05*(y_max_4-y_min_4), y_min_4, y_max_4];
drawAxes(ax2, X_axis_4, Y_axis_4, 0.04, 0, 3);

hold(ax2, 'off');

% Save figure
fig2_path = fullfile(output_dir, 'sc_spatial_selectivity_4locations.pdf');
pdfSave(fig2_path, fig2.Position(3:4)/72, fig2);
fprintf('Saved: %s\n', fig2_path);
fprintf('Y-axis limits: [%.2f, %.2f]\n', y_min_4, y_max_4);

close(fig2);

fprintf('\nDone.\n');

end

%% Helper Function: Wilson Score Confidence Interval
function [ci_lo, ci_hi] = wilson_score_ci(k, n, alpha)
% wilson_score_ci - Compute Wilson score confidence interval for a proportion.
%
% INPUTS:
%   k     - Number of successes
%   n     - Total number of trials
%   alpha - Significance level (default 0.05 for 95% CI)
%
% OUTPUTS:
%   ci_lo - Lower bound of CI
%   ci_hi - Upper bound of CI

if nargin < 3
    alpha = 0.05;
end

if n == 0
    ci_lo = 0;
    ci_hi = 0;
    return;
end

p_hat = k / n;
z = norminv(1 - alpha/2);  % z for 95% CI = 1.96

% Wilson score interval
denominator = 1 + z^2/n;
center = (p_hat + z^2/(2*n)) / denominator;
margin = (z / denominator) * sqrt(p_hat*(1-p_hat)/n + z^2/(4*n^2));

ci_lo = max(0, center - margin);
ci_hi = min(1, center + margin);

end
