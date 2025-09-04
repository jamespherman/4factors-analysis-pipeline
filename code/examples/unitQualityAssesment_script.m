
clear
close all
sc = load('Feynman_08_20_2024_SC.mat');
fileDate = '08_20_2024';

% 1: Unit Quality Assessment

min_isi_ms = 1.5;
min_firing_rate_hz = 1.0;
nClusters = sc.nClusters;
goodUnitIndices = [];

% loop thru each cluster to find good units
for i = 1:nClusters
    spike_times_s = sc.clusterTimes{i};
    
    % skip if cluster is empty
    if isempty(spike_times_s) || length(spike_times_s) < 2
        continue;
    end
    
    % calculate ISIs in ms
    isi = diff(spike_times_s) * 1000;
    
    % calculate refractory period violations
    refractory_violations = sum(isi < min_isi_ms);
    violation_rate_percent = (refractory_violations / length(isi)) * 100;
    
    % calculate mean firing rate
    session_duration_s = sc.spikeTimes(end); % use total recording duration
    mean_firing_rate = length(spike_times_s) / session_duration_s;
    
    % check if unit meets quality criteria
    if violation_rate_percent < 1.5 && mean_firing_rate > min_firing_rate_hz
        goodUnitIndices = [goodUnitIndices; i];
    end
end

% 2: Prepare Data for Plotting

% define bad unit indices
allUnitIndices = 1:nClusters;
badUnitIndices = setdiff(allUnitIndices, goodUnitIndices)';

% define cluster groups with specified colors (RGB values normalized to 1)
clusterGroups = { ...
    struct('indices', goodUnitIndices, 'color', [64, 176, 166]/255), ... % green
    struct('indices', badUnitIndices, 'color', [230, 97, 90]/255)  ... % red
};

% define event data and names for the PSTHs from the saccade task

sc.gSac = initializeGSacContrast(sc.trialInfo, sc.eventTimes);
eventTimes = {sc.gSac.targetOnTime(~isnan(sc.gSac.targetOnTime))};
eventNames = {'Target On'};

% 3: Generate Plot

makePSTHGrid(sc.clusterTimes, eventTimes, ...
    'ClusterGroups', clusterGroups, ...
    'EventNames', eventNames, ...
    'TimeRange', [-0.1, 0.5], ...
    'BinWidth', 0.025, ...
    'Title', 'PSTH for Good (Green) vs. Bad (Red) Units');