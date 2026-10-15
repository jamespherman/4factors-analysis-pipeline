function generate_sc_neuron_report()
% generate_sc_neuron_report - Creates a comprehensive report about SC neurons.
%
% This script generates a markdown report with:
%   1. Neuron counts by subject (Feynman, Newton) and type (classic, interneuron)
%   2. Significant target onset responses by target location (1-4) for each type
%   3. Significant modulation by each of the 4 factors at visual epoch
%
% The report is saved to docs/SC_Neuron_Report.md
%
% Author: Claude Code
% Date: 2026-02-01

%% Setup Paths
project_root = fullfile(findOneDrive, 'Code', '4factors-analysis-pipeline');
addpath(fullfile(project_root, 'code', 'utils'));

%% Load Aggregated Data
fprintf('Loading aggregated analysis data...\n');
agg_data_path = fullfile(project_root, 'data', 'processed', 'aggregated_analysis_data.mat');
if ~exist(agg_data_path, 'file')
    error('Aggregated data not found at: %s\nRun aggregate_analysis_results first.', agg_data_path);
end
load(agg_data_path, 'aggregated_sc_data', 'session_ids');

%% Extract Neuron Metrics Table
if ~isfield(aggregated_sc_data, 'neuron_metrics_table') || ...
        isempty(aggregated_sc_data.neuron_metrics_table)
    error('neuron_metrics_table not found in aggregated_sc_data.');
end
nmt = aggregated_sc_data.neuron_metrics_table;
fprintf('Loaded %d total SC neurons.\n', height(nmt));

%% Section 1: Count Neurons by Subject and Type
fprintf('Counting neurons by subject and type...\n');

subjects = {'Feynman', 'Newton'};
subject_counts = struct();

for i_subj = 1:length(subjects)
    subj = subjects{i_subj};

    % Find neurons from this subject (session_id starts with subject name)
    subj_mask = startsWith(nmt.session_id, subj);

    % Count by neuron type
    subject_counts.(subj).n_classic = sum(strcmp(nmt.neuron_type(subj_mask), 'classic'));
    subject_counts.(subj).n_interneuron = sum(strcmp(nmt.neuron_type(subj_mask), 'interneuron'));
    subject_counts.(subj).n_selected = subject_counts.(subj).n_classic + subject_counts.(subj).n_interneuron;

    % Count excluded (neuron_type starts with 'excluded')
    subject_counts.(subj).n_excluded = sum(startsWith(nmt.neuron_type(subj_mask), 'excluded'));
    subject_counts.(subj).n_total = sum(subj_mask);

    % Count sessions
    unique_sessions = unique(nmt.session_id(subj_mask));
    subject_counts.(subj).n_sessions = length(unique_sessions);

    fprintf('  %s: %d classic, %d interneuron, %d excluded (%d sessions)\n', ...
        subj, subject_counts.(subj).n_classic, subject_counts.(subj).n_interneuron, ...
        subject_counts.(subj).n_excluded, subject_counts.(subj).n_sessions);
end

% Totals
total_classic = 0;
total_interneuron = 0;
total_excluded = 0;
total_neurons = 0;
total_sessions = length(session_ids.sc);

for i_subj = 1:length(subjects)
    subj = subjects{i_subj};
    total_classic = total_classic + subject_counts.(subj).n_classic;
    total_interneuron = total_interneuron + subject_counts.(subj).n_interneuron;
    total_excluded = total_excluded + subject_counts.(subj).n_excluded;
    total_neurons = total_neurons + subject_counts.(subj).n_total;
end

%% Section 1b: Count Neurons by Session
fprintf('Counting neurons by session...\n');

% Get all unique session IDs and sort them
all_sessions = unique(nmt.session_id);
session_counts = struct();

for i_sess = 1:length(all_sessions)
    sess_id = all_sessions{i_sess};

    % Create valid field name (replace dashes/spaces with underscores)
    field_name = matlab.lang.makeValidName(sess_id);

    % Find neurons from this session
    sess_mask = strcmp(nmt.session_id, sess_id);

    % Extract subject name (first part before underscore)
    parts = strsplit(sess_id, '_');
    subject_name = parts{1};

    % Extract date (parts 2-4)
    if length(parts) >= 4
        date_str = sprintf('%s/%s/%s', parts{2}, parts{3}, parts{4});
    else
        date_str = 'unknown';
    end

    % Count by neuron type
    session_counts.(field_name).session_id = sess_id;
    session_counts.(field_name).subject = subject_name;
    session_counts.(field_name).date = date_str;
    session_counts.(field_name).n_classic = sum(strcmp(nmt.neuron_type(sess_mask), 'classic'));
    session_counts.(field_name).n_interneuron = sum(strcmp(nmt.neuron_type(sess_mask), 'interneuron'));
    session_counts.(field_name).n_selected = session_counts.(field_name).n_classic + session_counts.(field_name).n_interneuron;
    session_counts.(field_name).n_excluded = sum(startsWith(nmt.neuron_type(sess_mask), 'excluded'));
    session_counts.(field_name).n_total = sum(sess_mask);

    fprintf('  %s: %d classic, %d interneuron, %d excluded\n', ...
        sess_id, session_counts.(field_name).n_classic, ...
        session_counts.(field_name).n_interneuron, ...
        session_counts.(field_name).n_excluded);
end

%% Section 2: Target Onset Responses by Location
% Load per-session data to get cached_friedman_results
fprintf('Loading per-session data for target response analysis...\n');

% Initialize counters for each location (1-4) and neuron type
% We'll track activations and suppressions separately
location_response = struct();
for loc = 1:4
    location_response(loc).classic_activated = 0;
    location_response(loc).classic_suppressed = 0;
    location_response(loc).classic_total = 0;
    location_response(loc).interneuron_activated = 0;
    location_response(loc).interneuron_suppressed = 0;
    location_response(loc).interneuron_total = 0;
end

% Define standard 4factors target angles (in degrees)
% Location 1 = 45 deg, Location 2 = 135 deg, Location 3 = 225 deg, Location 4 = 315 deg
standard_angles = [45, 135, 225, 315];

% Process each SC session
for i_sess = 1:length(session_ids.sc)
    session_id = session_ids.sc{i_sess};

    session_data_path = fullfile(project_root, 'data', 'processed', ...
        session_id, [session_id '_session_data.mat']);

    if ~exist(session_data_path, 'file')
        fprintf('  Skipping %s (file not found)\n', session_id);
        continue;
    end

    data = load(session_data_path, 'session_data');
    if ~isfield(data.session_data, 'cached_friedman_results')
        fprintf('  Skipping %s (no cached_friedman_results)\n', session_id);
        continue;
    end

    friedman = data.session_data.cached_friedman_results;

    % Get neuron classifications for this session
    session_mask = strcmp(nmt.session_id, session_id);
    session_neuron_types = nmt.neuron_type(session_mask);

    % Validate dimensions
    n_clusters = size(friedman.sig_by_location, 1);
    if length(session_neuron_types) ~= n_clusters
        fprintf('  Warning: dimension mismatch for %s\n', session_id);
        continue;
    end

    % Map location_angles to standard locations (1-4)
    location_angles = friedman.location_angles;
    n_groups = length(location_angles);

    % Process each neuron
    for i_neuron = 1:n_clusters
        neuron_type = session_neuron_types{i_neuron};

        % Skip excluded neurons
        if startsWith(neuron_type, 'excluded')
            continue;
        end

        is_classic = strcmp(neuron_type, 'classic');
        is_interneuron = strcmp(neuron_type, 'interneuron');

        if ~is_classic && ~is_interneuron
            continue;
        end

        % Check each location group
        for i_group = 1:n_groups
            angle = location_angles(i_group);

            % Map angle to location index (1-4)
            [~, loc_idx] = min(abs(standard_angles - angle));

            % Only consider 4factors locations (not jph trials which may have different angles)
            if abs(standard_angles(loc_idx) - angle) > 5
                continue; % Skip non-standard angles
            end

            % Check for activation or suppression (visual epoch = index 1)
            % sig_by_location: [nClusters x nGroups x 3] for vis, delay, sac
            has_activation = friedman.sig_by_location(i_neuron, i_group, 1);
            has_suppression = false;
            if isfield(friedman, 'suppression_by_location')
                has_suppression = friedman.suppression_by_location(i_neuron, i_group, 1);
            end

            % Update counters
            if is_classic
                location_response(loc_idx).classic_total = location_response(loc_idx).classic_total + 1;
                if has_activation
                    location_response(loc_idx).classic_activated = location_response(loc_idx).classic_activated + 1;
                end
                if has_suppression
                    location_response(loc_idx).classic_suppressed = location_response(loc_idx).classic_suppressed + 1;
                end
            else % interneuron
                location_response(loc_idx).interneuron_total = location_response(loc_idx).interneuron_total + 1;
                if has_activation
                    location_response(loc_idx).interneuron_activated = location_response(loc_idx).interneuron_activated + 1;
                end
                if has_suppression
                    location_response(loc_idx).interneuron_suppressed = location_response(loc_idx).interneuron_suppressed + 1;
                end
            end
        end
    end
end

%% Section 3: 4-Factors Modulation (Visual Epoch)
fprintf('Analyzing 4-factors modulation at visual epoch...\n');

factors = {'reward', 'salience', 'identity', 'probability'};
factor_modulation = struct();

% Get masks for classic and interneuron
is_classic_mask = strcmp(nmt.neuron_type, 'classic');
is_interneuron_mask = strcmp(nmt.neuron_type, 'interneuron');

for i_factor = 1:length(factors)
    factor_name = factors{i_factor};
    p_col = sprintf('p_visual_%s', factor_name);

    if ~ismember(p_col, nmt.Properties.VariableNames)
        fprintf('  Warning: column %s not found\n', p_col);
        factor_modulation.(factor_name).classic_sig = 0;
        factor_modulation.(factor_name).classic_total = sum(is_classic_mask);
        factor_modulation.(factor_name).interneuron_sig = 0;
        factor_modulation.(factor_name).interneuron_total = sum(is_interneuron_mask);
        continue;
    end

    p_vals = nmt.(p_col);

    % Classic neurons
    classic_valid = is_classic_mask & ~isnan(p_vals);
    factor_modulation.(factor_name).classic_sig = sum(p_vals(classic_valid) < 0.05);
    factor_modulation.(factor_name).classic_total = sum(classic_valid);

    % Interneurons
    intern_valid = is_interneuron_mask & ~isnan(p_vals);
    factor_modulation.(factor_name).interneuron_sig = sum(p_vals(intern_valid) < 0.05);
    factor_modulation.(factor_name).interneuron_total = sum(intern_valid);

    fprintf('  %s: %d/%d classic sig, %d/%d interneuron sig\n', ...
        factor_name, ...
        factor_modulation.(factor_name).classic_sig, factor_modulation.(factor_name).classic_total, ...
        factor_modulation.(factor_name).interneuron_sig, factor_modulation.(factor_name).interneuron_total);
end

%% Generate Markdown Report
fprintf('Generating markdown report...\n');

report_path = fullfile(project_root, 'docs', 'SC_Neuron_Report.md');
fid = fopen(report_path, 'w');

% Header
fprintf(fid, '# SC Neuron Population Report\n\n');
fprintf(fid, '_Generated: %s_\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% Summary
fprintf(fid, '## Summary\n\n');
fprintf(fid, '- **Total SC neurons recorded**: %d\n', total_neurons);
fprintf(fid, '- **Total selected (classic + interneuron)**: %d\n', total_classic + total_interneuron);
fprintf(fid, '- **Classic SC neurons**: %d\n', total_classic);
fprintf(fid, '- **Putative interneurons**: %d\n', total_interneuron);
fprintf(fid, '- **Excluded**: %d\n', total_excluded);
fprintf(fid, '- **Total sessions**: %d\n\n', total_sessions);

% Neuron counts by subject
fprintf(fid, '## Neuron Counts by Subject\n\n');
fprintf(fid, '| Subject | Classic | Interneuron | Total Selected | Excluded | Total Recorded | Sessions |\n');
fprintf(fid, '|---------|---------|-------------|----------------|----------|----------------|----------|\n');

for i_subj = 1:length(subjects)
    subj = subjects{i_subj};
    sc = subject_counts.(subj);
    fprintf(fid, '| %s | %d | %d | %d | %d | %d | %d |\n', ...
        subj, sc.n_classic, sc.n_interneuron, sc.n_selected, ...
        sc.n_excluded, sc.n_total, sc.n_sessions);
end

% Add totals row
fprintf(fid, '| **Total** | **%d** | **%d** | **%d** | **%d** | **%d** | **%d** |\n\n', ...
    total_classic, total_interneuron, total_classic + total_interneuron, ...
    total_excluded, total_neurons, total_sessions);

% Neuron counts by session
fprintf(fid, '## Neuron Counts by Session\n\n');
fprintf(fid, '| Session | Subject | Date | Classic | Interneuron | Total Selected | Excluded | Total Recorded |\n');
fprintf(fid, '|---------|---------|------|---------|-------------|----------------|----------|----------------|\n');

session_field_names = fieldnames(session_counts);
for i_sess = 1:length(session_field_names)
    field_name = session_field_names{i_sess};
    sc = session_counts.(field_name);
    fprintf(fid, '| %s | %s | %s | %d | %d | %d | %d | %d |\n', ...
        sc.session_id, sc.subject, sc.date, ...
        sc.n_classic, sc.n_interneuron, sc.n_selected, ...
        sc.n_excluded, sc.n_total);
end
fprintf(fid, '\n');

% Target onset responses by location
fprintf(fid, '## Target Onset Responses by Location (Visual Epoch)\n\n');
fprintf(fid, 'Counts of neurons with significant activation or suppression at target onset, by target location.\n\n');
fprintf(fid, '| Location | Angle | Classic Activated | Classic Suppressed | Interneuron Activated | Interneuron Suppressed |\n');
fprintf(fid, '|----------|-------|-------------------|--------------------|-----------------------|------------------------|\n');

for loc = 1:4
    lr = location_response(loc);

    % Calculate percentages
    if lr.classic_total > 0
        classic_act_pct = 100 * lr.classic_activated / lr.classic_total;
        classic_sup_pct = 100 * lr.classic_suppressed / lr.classic_total;
    else
        classic_act_pct = 0;
        classic_sup_pct = 0;
    end

    if lr.interneuron_total > 0
        intern_act_pct = 100 * lr.interneuron_activated / lr.interneuron_total;
        intern_sup_pct = 100 * lr.interneuron_suppressed / lr.interneuron_total;
    else
        intern_act_pct = 0;
        intern_sup_pct = 0;
    end

    fprintf(fid, '| %d | %d° | %d (%.1f%%) | %d (%.1f%%) | %d (%.1f%%) | %d (%.1f%%) |\n', ...
        loc, standard_angles(loc), ...
        lr.classic_activated, classic_act_pct, ...
        lr.classic_suppressed, classic_sup_pct, ...
        lr.interneuron_activated, intern_act_pct, ...
        lr.interneuron_suppressed, intern_sup_pct);
end

fprintf(fid, '\n_Note: Each neuron is counted once per location where it was tested. ');
fprintf(fid, 'Percentages are relative to the number of neurons tested at each location._\n\n');

% 4-Factors modulation
fprintf(fid, '## 4-Factors Modulation (Visual Epoch)\n\n');
fprintf(fid, 'Proportion of neurons with significant (p < 0.05) selectivity for each factor during the visual response epoch.\n\n');
fprintf(fid, '| Factor | Classic Sig | Classic Total | Classic %% | Interneuron Sig | Interneuron Total | Interneuron %% |\n');
fprintf(fid, '|--------|-------------|---------------|-----------|-----------------|-------------------|---------------|\n');

for i_factor = 1:length(factors)
    factor_name = factors{i_factor};
    fm = factor_modulation.(factor_name);

    if fm.classic_total > 0
        classic_pct = 100 * fm.classic_sig / fm.classic_total;
    else
        classic_pct = 0;
    end

    if fm.interneuron_total > 0
        intern_pct = 100 * fm.interneuron_sig / fm.interneuron_total;
    else
        intern_pct = 0;
    end

    % Capitalize factor name for display
    factor_display = [upper(factor_name(1)), factor_name(2:end)];

    fprintf(fid, '| %s | %d | %d | %.1f%% | %d | %d | %.1f%% |\n', ...
        factor_display, fm.classic_sig, fm.classic_total, classic_pct, ...
        fm.interneuron_sig, fm.interneuron_total, intern_pct);
end

fprintf(fid, '\n');

% Methods note
fprintf(fid, '## Methods\n\n');
fprintf(fid, '### Neuron Classification\n\n');
fprintf(fid, '**Classic SC neurons** meet all criteria:\n');
fprintf(fid, '- Significant increase in visual/delay/saccade epoch vs baseline for ≥1 contralateral location\n');
fprintf(fid, '- Mean FR > 3.0 sp/s\n');
fprintf(fid, '- NOT significantly activated at opposing locations (180° apart)\n');
fprintf(fid, '- Not sparse (≤80%% of 100ms bins empty)\n\n');

fprintf(fid, '**Putative interneurons** meet:\n');
fprintf(fid, '- Mean FR > 3.0 sp/s\n');
fprintf(fid, '- Task-modulated (activation or suppression somewhere)\n');
fprintf(fid, '- Either opposing activation OR contralateral suppression\n\n');

fprintf(fid, '**Excluded** neurons fail to meet basic criteria (low FR, sparse, or not task-modulated).\n\n');

fprintf(fid, '### Statistical Tests\n\n');
fprintf(fid, '- Target onset responses: Friedman test with post-hoc Wilcoxon signed-rank tests (Bonferroni corrected)\n');
fprintf(fid, '- 4-Factors modulation: Window-based ROC analysis comparing high vs low conditions\n');

fclose(fid);

fprintf('Report saved to: %s\n', report_path);
fprintf('Done.\n');

end
