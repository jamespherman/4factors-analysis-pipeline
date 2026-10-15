function generate_snc_neuron_report()
% generate_snc_neuron_report - Creates a comprehensive report about SNc neurons.
%
% This script generates a markdown report with:
%   1. Neuron counts by subject and by session
%   2. Selected (modulated) vs excluded counts
%   3. Breakdown of selected neurons by baseline firing rate (<20 sp/s threshold)
%
% The report is saved to docs/SNc_Neuron_Report.md
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
load(agg_data_path, 'aggregated_snc_data', 'session_ids');

%% Extract Neuron Metrics Table
if ~isfield(aggregated_snc_data, 'neuron_metrics_table') || ...
        isempty(aggregated_snc_data.neuron_metrics_table)
    error('neuron_metrics_table not found in aggregated_snc_data.');
end
nmt = aggregated_snc_data.neuron_metrics_table;
fprintf('Loaded %d total SNc neurons.\n', height(nmt));

%% Section 1: Count Neurons by Subject
fprintf('Counting neurons by subject...\n');

% Get unique subjects from session IDs
all_session_ids = unique(nmt.session_id);
subjects_from_sessions = cellfun(@(x) strsplit(x, '_'), all_session_ids, 'UniformOutput', false);
subjects_from_sessions = cellfun(@(x) x{1}, subjects_from_sessions, 'UniformOutput', false);
subjects = unique(subjects_from_sessions);

subject_counts = struct();
total_selected = 0;
total_excluded = 0;
total_neurons = 0;

for i_subj = 1:length(subjects)
    subj = subjects{i_subj};

    % Find neurons from this subject
    subj_mask = startsWith(nmt.session_id, subj);

    % Count by neuron type (modulated vs excluded for SNc)
    subject_counts.(subj).n_modulated = sum(strcmp(nmt.neuron_type(subj_mask), 'modulated'));
    subject_counts.(subj).n_excluded = sum(startsWith(nmt.neuron_type(subj_mask), 'excluded'));
    subject_counts.(subj).n_total = sum(subj_mask);

    % Count sessions
    unique_sessions = unique(nmt.session_id(subj_mask));
    subject_counts.(subj).n_sessions = length(unique_sessions);

    % For selected neurons, count those with baseline_fr < 20 sp/s
    selected_mask = subj_mask & strcmp(nmt.neuron_type, 'modulated');
    baseline_frs = nmt.baseline_fr(selected_mask);
    subject_counts.(subj).n_low_baseline = sum(baseline_frs < 20);
    subject_counts.(subj).n_high_baseline = sum(baseline_frs >= 20);

    % Update totals
    total_selected = total_selected + subject_counts.(subj).n_modulated;
    total_excluded = total_excluded + subject_counts.(subj).n_excluded;
    total_neurons = total_neurons + subject_counts.(subj).n_total;

    fprintf('  %s: %d modulated, %d excluded (%d sessions)\n', ...
        subj, subject_counts.(subj).n_modulated, ...
        subject_counts.(subj).n_excluded, subject_counts.(subj).n_sessions);
end

total_sessions = length(session_ids.snc);

%% Section 2: Count Neurons by Session
fprintf('Counting neurons by session...\n');

session_counts = struct();

for i_sess = 1:length(all_session_ids)
    sess_id = all_session_ids{i_sess};
    field_name = matlab.lang.makeValidName(sess_id);

    % Find neurons from this session
    sess_mask = strcmp(nmt.session_id, sess_id);

    % Extract subject and date
    parts = strsplit(sess_id, '_');
    subject_name = parts{1};
    if length(parts) >= 4
        date_str = sprintf('%s/%s/%s', parts{2}, parts{3}, parts{4});
    else
        date_str = 'unknown';
    end

    % Count by type
    session_counts.(field_name).session_id = sess_id;
    session_counts.(field_name).subject = subject_name;
    session_counts.(field_name).date = date_str;
    session_counts.(field_name).n_modulated = sum(strcmp(nmt.neuron_type(sess_mask), 'modulated'));
    session_counts.(field_name).n_excluded = sum(startsWith(nmt.neuron_type(sess_mask), 'excluded'));
    session_counts.(field_name).n_total = sum(sess_mask);

    % For selected neurons, count those with baseline_fr < 20 sp/s
    selected_mask = sess_mask & strcmp(nmt.neuron_type, 'modulated');
    baseline_frs = nmt.baseline_fr(selected_mask);
    session_counts.(field_name).n_low_baseline = sum(baseline_frs < 20);
    session_counts.(field_name).n_high_baseline = sum(baseline_frs >= 20);

    fprintf('  %s: %d modulated (%d low FR, %d high FR), %d excluded\n', ...
        sess_id, session_counts.(field_name).n_modulated, ...
        session_counts.(field_name).n_low_baseline, ...
        session_counts.(field_name).n_high_baseline, ...
        session_counts.(field_name).n_excluded);
end

%% Section 3: Overall baseline FR statistics for selected neurons
fprintf('Computing baseline FR statistics for selected neurons...\n');

selected_mask = strcmp(nmt.neuron_type, 'modulated');
selected_baseline_frs = nmt.baseline_fr(selected_mask);

% Remove NaNs for statistics
valid_baseline = selected_baseline_frs(~isnan(selected_baseline_frs));

n_below_5 = sum(valid_baseline < 5);
n_5_to_10 = sum(valid_baseline >= 5 & valid_baseline < 10);
n_10_to_20 = sum(valid_baseline >= 10 & valid_baseline < 20);
n_20_plus = sum(valid_baseline >= 20);

fprintf('  < 5 sp/s: %d\n', n_below_5);
fprintf('  5-10 sp/s: %d\n', n_5_to_10);
fprintf('  10-20 sp/s: %d\n', n_10_to_20);
fprintf('  >= 20 sp/s: %d\n', n_20_plus);

%% Generate Markdown Report
fprintf('Generating markdown report...\n');

report_path = fullfile(project_root, 'docs', 'SNc_Neuron_Report.md');
fid = fopen(report_path, 'w');

% Header
fprintf(fid, '# SNc Neuron Population Report\n\n');
fprintf(fid, '_Generated: %s_\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

% Summary
fprintf(fid, '## Summary\n\n');
fprintf(fid, '- **Total SNc neurons recorded**: %d\n', total_neurons);
fprintf(fid, '- **Selected (modulated)**: %d (%.1f%%)\n', total_selected, 100*total_selected/total_neurons);
fprintf(fid, '- **Excluded**: %d (%.1f%%)\n', total_excluded, 100*total_excluded/total_neurons);
fprintf(fid, '- **Total sessions**: %d\n\n', total_sessions);

% Neuron counts by subject
fprintf(fid, '## Neuron Counts by Subject\n\n');
fprintf(fid, '| Subject | Modulated | Excluded | Total | Sessions |\n');
fprintf(fid, '|---------|-----------|----------|-------|----------|\n');

for i_subj = 1:length(subjects)
    subj = subjects{i_subj};
    sc = subject_counts.(subj);
    fprintf(fid, '| %s | %d | %d | %d | %d |\n', ...
        subj, sc.n_modulated, sc.n_excluded, sc.n_total, sc.n_sessions);
end

if length(subjects) > 1
    fprintf(fid, '| **Total** | **%d** | **%d** | **%d** | **%d** |\n\n', ...
        total_selected, total_excluded, total_neurons, total_sessions);
else
    fprintf(fid, '\n');
end

% Neuron counts by session
fprintf(fid, '## Neuron Counts by Session\n\n');
fprintf(fid, '| Session | Subject | Date | Modulated | Excluded | Total |\n');
fprintf(fid, '|---------|---------|------|-----------|----------|-------|\n');

session_field_names = fieldnames(session_counts);
for i_sess = 1:length(session_field_names)
    field_name = session_field_names{i_sess};
    sc = session_counts.(field_name);
    fprintf(fid, '| %s | %s | %s | %d | %d | %d |\n', ...
        sc.session_id, sc.subject, sc.date, ...
        sc.n_modulated, sc.n_excluded, sc.n_total);
end
fprintf(fid, '\n');

% Baseline FR breakdown for selected neurons
fprintf(fid, '## Baseline Firing Rate Distribution (Selected Neurons)\n\n');
fprintf(fid, 'Distribution of baseline firing rates among the %d selected (modulated) SNc neurons.\n\n', total_selected);
fprintf(fid, '| Baseline FR Range | Count | Percentage |\n');
fprintf(fid, '|-------------------|-------|------------|\n');
fprintf(fid, '| < 5 sp/s | %d | %.1f%% |\n', n_below_5, 100*n_below_5/length(valid_baseline));
fprintf(fid, '| 5-10 sp/s | %d | %.1f%% |\n', n_5_to_10, 100*n_5_to_10/length(valid_baseline));
fprintf(fid, '| 10-20 sp/s | %d | %.1f%% |\n', n_10_to_20, 100*n_10_to_20/length(valid_baseline));
fprintf(fid, '| >= 20 sp/s | %d | %.1f%% |\n', n_20_plus, 100*n_20_plus/length(valid_baseline));
fprintf(fid, '| **< 20 sp/s (Total)** | **%d** | **%.1f%%** |\n\n', ...
    n_below_5 + n_5_to_10 + n_10_to_20, ...
    100*(n_below_5 + n_5_to_10 + n_10_to_20)/length(valid_baseline));

% Per-session breakdown with baseline FR
fprintf(fid, '## Baseline Firing Rate by Session (Selected Neurons)\n\n');
fprintf(fid, '| Session | Modulated | < 20 sp/s | >= 20 sp/s | %% Low FR |\n');
fprintf(fid, '|---------|-----------|-----------|------------|----------|\n');

for i_sess = 1:length(session_field_names)
    field_name = session_field_names{i_sess};
    sc = session_counts.(field_name);
    if sc.n_modulated > 0
        pct_low = 100 * sc.n_low_baseline / sc.n_modulated;
    else
        pct_low = 0;
    end
    fprintf(fid, '| %s | %d | %d | %d | %.1f%% |\n', ...
        sc.session_id, sc.n_modulated, sc.n_low_baseline, sc.n_high_baseline, pct_low);
end
fprintf(fid, '\n');

% Summary statistics
fprintf(fid, '## Baseline FR Summary Statistics (Selected Neurons)\n\n');
fprintf(fid, '| Statistic | Value |\n');
fprintf(fid, '|-----------|-------|\n');
fprintf(fid, '| Mean | %.2f sp/s |\n', mean(valid_baseline));
fprintf(fid, '| Median | %.2f sp/s |\n', median(valid_baseline));
fprintf(fid, '| Std Dev | %.2f sp/s |\n', std(valid_baseline));
fprintf(fid, '| Min | %.2f sp/s |\n', min(valid_baseline));
fprintf(fid, '| Max | %.2f sp/s |\n', max(valid_baseline));
fprintf(fid, '\n');

% Methods note
fprintf(fid, '## Methods\n\n');
fprintf(fid, '### SNc Neuron Classification\n\n');
fprintf(fid, '**Modulated SNc neurons** meet all criteria:\n');
fprintf(fid, '- Significant change in activity (increase or decrease) for at least 3 consecutive samples\n');
fprintf(fid, '- Modulation detected at any of: fixation onset, target onset, fixation offset, peri-saccade, or reward\n');
fprintf(fid, '- Mean FR > 1 sp/s (computed over active period)\n');
fprintf(fid, '- Not sparse (<=80%% of 100ms bins empty)\n\n');

fprintf(fid, '**Excluded** neurons fail to meet one or more of the above criteria.\n\n');

fprintf(fid, '### Baseline Firing Rate\n\n');
fprintf(fid, 'Baseline FR is computed during the pre-target baseline epoch (-75 to +25 ms relative to target onset).\n');
fprintf(fid, 'The 20 sp/s threshold is used to identify putative dopamine neurons, as dopamine neurons typically\n');
fprintf(fid, 'have low baseline firing rates (<10 sp/s), while GABAergic interneurons often have higher rates.\n');

fclose(fid);

fprintf('Report saved to: %s\n', report_path);
fprintf('Done.\n');

end
