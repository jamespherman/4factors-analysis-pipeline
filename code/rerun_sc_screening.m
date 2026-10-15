% rerun_sc_screening.m - Re-run SC neuron screening with updated classification criteria
% This script clears cached Friedman results and re-runs screening for all SC sessions
% to apply the new interneuron detection criteria:
%   - Opposing activation check uses only Visual and Saccade epochs (not Delay)
%   - Suppression-based interneuron detection added
%
% Updates both session_data (cached Friedman results) and analysis_results (screening_info)

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_root, 'code'));
addpath(fullfile(project_root, 'code', 'utils'));

% Find all SC sessions dynamically
processed_dir = fullfile(project_root, 'data', 'processed');
all_dirs = dir(processed_dir);
sc_sessions = {};
for i = 1:length(all_dirs)
    if all_dirs(i).isdir && endsWith(all_dirs(i).name, '_SC')
        sc_sessions{end+1} = all_dirs(i).name;
    end
end
sc_sessions = sc_sessions';
fprintf('Found %d SC sessions.\n', length(sc_sessions));

% Load condition definitions
condition_defs = define_task_conditions();

fprintf('=== Re-running SC Neuron Screening with Updated Criteria ===\n');
fprintf('New criteria:\n');
fprintf('  - Opposing activation: Visual + Saccade only (not Delay)\n');
fprintf('  - Suppression-based interneuron detection enabled\n\n');

% Track results
results = struct();
results.session_name = {};
results.n_classic = [];
results.n_interneuron = [];
results.n_excluded = [];
results.n_total = [];

for i_sess = 1:length(sc_sessions)
    session_name = sc_sessions{i_sess};
    fprintf('\n--- Processing %s (%d/%d) ---\n', session_name, i_sess, length(sc_sessions));

    % Load session data
    session_file = fullfile(project_root, 'data', 'processed', session_name, ...
        [session_name, '_session_data.mat']);

    if ~exist(session_file, 'file')
        fprintf('WARNING: Session file not found: %s\n', session_file);
        continue;
    end

    loaded = load(session_file);
    session_data = loaded.session_data;

    % Clear cached Friedman results to force recomputation
    if isfield(session_data, 'cached_friedman_results')
        fprintf('Clearing cached Friedman results...\n');
        session_data = rmfield(session_data, 'cached_friedman_results');
    end

    % Run screening
    [selected_neurons, sig_epoch_comparison, scSide, neuron_class, session_data] = ...
        screen_sc_neurons(session_data, project_root, condition_defs, session_name);

    % Count results
    n_classic = sum(strcmp(neuron_class, 'classic'));
    n_interneuron = sum(strcmp(neuron_class, 'interneuron'));
    n_excluded = sum(strcmp(neuron_class, 'excluded'));
    n_total = length(neuron_class);

    % Store results
    results.session_name{end+1} = session_name;
    results.n_classic(end+1) = n_classic;
    results.n_interneuron(end+1) = n_interneuron;
    results.n_excluded(end+1) = n_excluded;
    results.n_total(end+1) = n_total;

    % Also update session_data.analysis with new classifications
    if isfield(session_data, 'analysis')
        % Update screening_info.neuron_class
        if isfield(session_data.analysis, 'screening_info')
            fprintf('Updating session_data.analysis.screening_info with new neuron classifications...\n');
            for i = 1:length(session_data.analysis.screening_info)
                if i <= length(neuron_class)
                    session_data.analysis.screening_info(i).neuron_class = neuron_class{i};
                    % Also update 'included' field based on neuron_class
                    if strcmp(neuron_class{i}, 'classic') || strcmp(neuron_class{i}, 'interneuron')
                        session_data.analysis.screening_info(i).included = true;
                    else
                        session_data.analysis.screening_info(i).included = false;
                    end
                end
            end
        end

        % Update selected_neurons array
        session_data.analysis.selected_neurons = selected_neurons;
        fprintf('Updated session_data.analysis.selected_neurons.\n');
    else
        fprintf('WARNING: No session_data.analysis found, skipping update.\n');
    end

    % Save updated session data with new cached results and updated screening_info
    save(session_file, 'session_data');
    fprintf('Saved updated session data.\n');
end

% Print summary
fprintf('\n\n========================================\n');
fprintf('        SC NEURON SCREENING SUMMARY\n');
fprintf('========================================\n\n');
fprintf('%-30s %8s %12s %10s %8s\n', 'Session', 'Classic', 'Interneuron', 'Excluded', 'Total');
fprintf('%s\n', repmat('-', 1, 72));

for i = 1:length(results.session_name)
    fprintf('%-30s %8d %12d %10d %8d\n', ...
        results.session_name{i}, ...
        results.n_classic(i), ...
        results.n_interneuron(i), ...
        results.n_excluded(i), ...
        results.n_total(i));
end

fprintf('%s\n', repmat('-', 1, 72));
fprintf('%-30s %8d %12d %10d %8d\n', ...
    'TOTAL', ...
    sum(results.n_classic), ...
    sum(results.n_interneuron), ...
    sum(results.n_excluded), ...
    sum(results.n_total));
fprintf('\n');

% Save results summary
summary_file = fullfile(project_root, 'data', 'sc_screening_summary.mat');
save(summary_file, 'results');
fprintf('Summary saved to: %s\n', summary_file);
