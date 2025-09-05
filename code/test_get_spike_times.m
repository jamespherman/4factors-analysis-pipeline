%% Test script for the data_handling.get_spike_times function
%
% This script tests the functionality of the get_spike_times function to
% ensure it correctly extracts spike times for given neuron IDs.

%% Setup
% Add the code directory to the MATLAB path.
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(script_dir);
fprintf('Setup complete. Paths configured.\n');

%% Test 1: Basic functionality
try
    % Create mock session_data
    session_data.spikes.times = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
    session_data.spikes.clusters = [1, 2, 1, 3, 2, 1];
    
    neuron_ids = [1, 2];
    
    % Call the function
    spike_times_cell = data_handling.get_spike_times(session_data, neuron_ids);
    
    % --- Verification ---
    assert(iscell(spike_times_cell), 'Output is not a cell array.');
    assert(numel(spike_times_cell) == 2, 'Output cell array has incorrect size.');
    assert(all(spike_times_cell{1} == [0.1, 0.3, 0.6]), 'Spike times for neuron 1 are incorrect.');
    assert(all(spike_times_cell{2} == [0.2, 0.5]), 'Spike times for neuron 2 are incorrect.');
    
    fprintf('[PASS] Basic functionality test passed.\n');
catch ME
    fprintf('[FAIL] Basic functionality test failed: %s\n', ME.message);
end

%% Test 2: Edge case - empty neuron_ids
try
    % Create mock session_data
    session_data.spikes.times = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
    session_data.spikes.clusters = [1, 2, 1, 3, 2, 1];
    
    neuron_ids = [];
    
    % Call the function
    spike_times_cell = data_handling.get_spike_times(session_data, neuron_ids);
    
    % --- Verification ---
    assert(iscell(spike_times_cell), 'Output is not a cell array for empty input.');
    assert(isempty(spike_times_cell), 'Output is not empty for empty input.');
    
    fprintf('[PASS] Edge case with empty neuron_ids passed.\n');
catch ME
    fprintf('[FAIL] Edge case with empty neuron_ids failed: %s\n', ME.message);
end

%% Test 3: Edge case - neuron_id not found
try
    % Create mock session_data
    session_data.spikes.times = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
    session_data.spikes.clusters = [1, 2, 1, 3, 2, 1];
    
    neuron_ids = [4, 1]; % Neuron 4 does not exist
    
    % Call the function
    spike_times_cell = data_handling.get_spike_times(session_data, neuron_ids);
    
    % --- Verification ---
    assert(iscell(spike_times_cell), 'Output is not a cell array.');
    assert(numel(spike_times_cell) == 2, 'Output cell array has incorrect size.');
    assert(isempty(spike_times_cell{1}), 'Cell for non-existent neuron is not empty.');
    assert(all(spike_times_cell{2} == [0.1, 0.3, 0.6]), 'Spike times for neuron 1 are incorrect.');
    
    fprintf('[PASS] Edge case with a non-existent neuron ID passed.\n');
catch ME
    fprintf('[FAIL] Edge case with a non-existent neuron ID failed: %s\n', ME.message);
end
