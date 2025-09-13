%% test_neuron_diagnostics.m
% This script serves as a step-by-step test of the neuron screening and
% diagnostic PDF workflow. It runs the necessary functions for a single,
% hardcoded session.

% Start timer
tic;

% In-line function to report timing
giveFeed = @(x)disp([num2str(toc) 's - ' x]);

%% Setup Paths
% Add the 'utils' directory to the path so that helper functions can be
% found.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

%% Setup
unique_id = 'Feynman_08_15_2025_SC'; % Hardcoded session
giveFeed(sprintf('Testing diagnostic workflow for session: %s', ...
    unique_id));

%% Step 1: Load session data
giveFeed('Step 1: Loading session data...');
one_drive_path = findOneDrive;
session_data_path = fullfile(one_drive_path, ...
    'Neuronal Data Analysis', unique_id, ...
    [unique_id '_session_data.mat']);

% Note: The following line will error if the data file does not exist.
% This is expected in a test environment without the actual data.
try
    data = load(session_data_path);
    session_data = data.session_data;
    giveFeed('Session data loaded successfully.');
catch ME
    warning('Could not load session data. Using dummy data instead.');
    warning(ME.message);
    % Create a dummy session_data struct for testing purposes
    session_data = struct('spikes', rand(100,1), 'neuron_ids', 1:5);
end

%% Step 2: Calculate and Store Diagnostic Metrics
giveFeed('Step 2: Calculating and storing diagnostic metrics...');

% Calculate baseline firing rates and add to session_data
session_data.metrics.baseline_frs = calculate_baseline_fr(session_data);

% Calculate waveform metrics and add to session_data
nClusters = height(session_data.spikes.cluster_info.cluster_id);
for i_cluster = 1:nClusters

    % get current unit's multi-channel waveform:
    mean_waveform = session_data.spikes.wfMeans{i_cluster};

    % find channel with max variance:
    [~,max_var_chan] = max(var(mean_waveform,[],2));

    % Note: Assuming a sampling rate of 30000 Hz
    session_data.metrics.wf_metrics(i_cluster, 1) = ...
        calculate_waveform_metrics(mean_waveform(max_var_chan,:), 30000);
end
giveFeed('Diagnostic metrics calculated and stored.');


%% Step 3: Run neuron screening
giveFeed('Step 3: Running neuron screening...');
selected_neurons = []; % Initialize empty

if contains(unique_id, 'SNc')
    giveFeed('Session is SNc type. Running screen_da_neurons...');
    % Assumes screen_da_neurons is in the path
    selected_neurons = screen_da_neurons(session_data, unique_id);
    giveFeed('DA neuron screening complete.');
elseif contains(unique_id, 'SC')
    giveFeed('Session is SC type. Running screen_sc_neurons...');

    % Store the results in the session_data structure
    if ~isfield(session_data, 'metadata')
        session_data.metadata = struct();
    end
    session_data.metadata.unique_id = unique_id;

    % screen_sc_neurons will now determine scSide and calculate 
    % significance
    [selected_neurons, sig_epoch_comp, scSide] = screen_sc_neurons( ...
        session_data);

    
    session_data.metadata.scSide = scSide;
    session_data.metadata.sig_epoch_comparison = sig_epoch_comp;

    giveFeed(sprintf(['SC neuron screening complete. Determined ' ...
        '        scSide: %s.'], scSide));
else
    error(['Unknown session type in unique_id. Cannot determine which ' ...
        '        screening function to run.']);
end

%% Step 4: Generate diagnostic PDF
giveFeed('Step 4: Generating diagnostic PDF...');
try
    generate_neuron_summary_pdf(session_data, selected_neurons, unique_id);
    giveFeed('Diagnostic PDF generation step complete.');
catch ME_pdf
    giveFeed('Error during PDF generation step.');
    rethrow(ME_pdf);
end

%% Done
giveFeed('Test complete. Check the generated PDF for session diagnostics.');
