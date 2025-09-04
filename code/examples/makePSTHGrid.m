function figHandle = makePSTHGrid(clusterTimes, eventTimes, varargin)
% psthGrid Generates PSTHs for specified groups of clusters.
%
% This function creates a tiled layout of peri-stimulus time histograms (PSTHs)
% for different groups of neuronal clusters. Each group can be assigned a
% specific color. Within each cluster's subplot, PSTHs for multiple event
% types are overlaid, distinguished by line style.
%
% Usage:
%   figHandle = makePSTHGrid(clusterTimes, eventTimes, 'ParameterName', ParameterValue, ...)
%
% Inputs:
%   clusterTimes      - Cell array of spike times for each cluster.
%   eventTimes        - Cell array of event time arrays, one cell per event type.
%
% Optional Name-Value Pair Arguments:
%   'ClusterGroups'   - Cell array of structs, where each struct defines a
%                       group with fields:
%                         .indices - Vector of cluster indices for the group.
%                         .color   - RGB triplet for the group's plot color.
%   'BinWidth'        - Width of PSTH bins in seconds (default: 0.025).
%   'TimeRange'       - [start, end] for the full PSTH time axis (default: [-0.1, 2]).
%   'EventNames'      - Cell array of strings for event names (for legend).
%   'Title'           - String for the main figure title.
%
% Outputs:
%   figHandle         - Handle to the figure containing the PSTH subplots.

%% parse input arguments
p = inputParser;
addRequired(p, 'clusterTimes', @iscell);
addRequired(p, 'eventTimes', @iscell);
addParameter(p, 'ClusterGroups', {}, @iscell);
addParameter(p, 'BinWidth', 0.025, @isscalar);
addParameter(p, 'TimeRange', [-0.1, 2], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'EventNames', {}, @iscell);
addParameter(p, 'Title', 'Cluster PSTHs');
parse(p, clusterTimes, eventTimes, varargin{:});

clusterGroups = p.Results.ClusterGroups;
binWidth = p.Results.BinWidth;
timeRange = p.Results.TimeRange;
eventNames = p.Results.EventNames;
figureTitle = p.Results.Title;

%% prepare data for plotting
% create a flat list of all clusters to plot and a map of their colors
allClusters = [];
clusterColorMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
for g = 1:length(clusterGroups)
    group = clusterGroups{g};
    allClusters = [1:max(max(clusterGroups{1}.indices), max(clusterGroups{2}.indices))];
    for k = 1:length(group.indices)
        clusterColorMap(group.indices(k)) = group.color;
    end
end

% create default event names if not provided
if isempty(eventNames)
    eventNames = cellstr(strcat('Event ', num2str((1:length(eventTimes))')));
end

%% perform calculations
% initialize variables for parallel processing
nEventTypes = length(eventTimes);
nTotalClusters = length(allClusters);
results = cell(nTotalClusters, nEventTypes);

% start parallel pool
if isempty(gcp('nocreate'))
    parpool('local');
end

% loop through all specified clusters and event types to calculate PSTHs
parfor i = 1:nTotalClusters
    clusterIndex = allClusters(i);
    spikeTimes = clusterTimes{clusterIndex};
    
    % create a temporary cell for the inner loop results
    tempResults = cell(1, nEventTypes);
    
    for j = 1:nEventTypes
        currentEventTimes = eventTimes{j};
        
        % align spikes to events and bin them
        [~, binnedCounts] = alignAndBinSpikes(spikeTimes, currentEventTimes, timeRange(1), timeRange(2), binWidth);
        
        % calculate confidence interval for the mean binned rate
        meanBinnedRate = mean(binnedCounts) / binWidth;
        meanBinnedRateCi = bootci(500, @(x) mean(x) / binWidth, binnedCounts);
        
        % store results for this event type
        tempResults{j} = struct('index', clusterIndex, ...
                                'meanBinnedRate', meanBinnedRate, 'meanBinnedRateCi', meanBinnedRateCi);
    end
    % assign the results for the cluster
    results(i, :) = tempResults;
end

%% generate plots
% determine the layout of the figure
nRows = ceil(sqrt(nTotalClusters));
nCols = ceil(nTotalClusters / nRows);

% create the figure and tiled layout
figHandle = figure('Name', figureTitle, 'Position', [100, 100, 300*nCols, 250*nRows]);
t = tiledlayout(figHandle, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

% define line styles to differentiate event types
lineStyles = {'-', '--', ':', '-.'};
alpha = 0.4; % transparency for confidence intervals

% iterate through clusters to create subplots
for i = 1:nTotalClusters
    clusterIndex = allClusters(i);
    groupColor = clusterColorMap(clusterIndex);
    
    % create a new tile for this cluster's PSTH
    ax = nexttile(t);
    hold on;
    
    % prepare handles for legend
    plotHandles = gobjects(1, nEventTypes);
    
    % plot PSTH with CI for all event types
    binCenters = timeRange(1)+binWidth/2 : binWidth : timeRange(2)-binWidth/2;
    for j = 1:nEventTypes
        res = results{i,j};
        ci = res.meanBinnedRateCi;
        currentLineStyle = lineStyles{mod(j-1, length(lineStyles)) + 1};
        
        % plot confidence interval shading
        x_fill = [binCenters, fliplr(binCenters)];
        y_fill = [ci(1,:), fliplr(ci(2,:))];
        fill(x_fill, y_fill, groupColor, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        
        % plot mean line with a unique line style
        plotHandles(j) = plot(binCenters, res.meanBinnedRate, 'Color', groupColor, 'LineWidth', 1.5, 'LineStyle', currentLineStyle);
    end
    hold off;
    
    % customize plot
    title(sprintf('Cluster %d', clusterIndex));
    xlabel('Time (s)');
    ylabel('Firing Rate (Hz)');
    xlim(timeRange);
    xline(0, ':', 'LineWidth', 1, 'Color', 'k');
    
    % add legend for the different event types
    legend(plotHandles, eventNames, 'Location', 'best', 'FontSize', 8, 'Box', 'off');
    
    % set axes properties
    set(gca, 'TickDir', 'Out', 'LineWidth', 1, 'FontSize', 10);
end

%% finalize figure
title(t, figureTitle);
set(figHandle, 'Color', 'w');

% add date annotation if fileDate exists
if evalin('base', 'exist(''fileDate'', ''var'')')
    filedateStr = evalin('base', 'fileDate');
    dateStr = strrep(filedateStr, '_', '-');
    annotation('textbox', [0.9, 0.01, 0.1, 0.05], 'String', dateStr, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'FontSize', 10);
end

end