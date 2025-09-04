function fig = plotPopulationPreferenceScOnly(aggSC, saveDir)
% PLOTPOPULATIONPREFERENCESCONLY Visualizes SC neuron population preference
%
%   fig = plotPopulationPreferenceScOnly(aggSC, saveDir)
%
% This function creates a visualization showing the proportion of SC neurons 
% preferring each condition across different task events and comparison types.
%
% Inputs:
%   aggSC    - Aggregated SC data structure from aggregateScSessions
%   saveDir  - Directory to save the output figure (optional)
%
% Output:
%   fig      - Handle to the created figure

% Handle optional input
if nargin < 2
    saveDir = pwd;
end

% Ensure the save directory exists
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% In-line function for timing feedback
tic;
giveFeed = @(x) disp([num2str(toc) ' - ' x]);

% Define color palette
palette = richColors;
palette = palette([6, 12], :);

% Get all events
events = fieldnames(aggSC);

% Group events by task type and count them
sacEvents = events(startsWith(events, 'sac_'));
attnEvents = events(startsWith(events, 'attn_'));
nSacEvents = length(sacEvents);
nAttnEvents = length(attnEvents);

% Dynamically determine comparisons for each task type
allComps = struct();
allComps.sac = {};
allComps.attn = {};

% Identify all unique comparisons for saccade task
if ~isempty(sacEvents)
    firstSacEvent = sacEvents{1};
    allCompsSac = fieldnames(aggSC.(firstSacEvent));
    allComps.sac = allCompsSac(startsWith(allCompsSac, 'raw_') & ...
        ~strcmp(allCompsSac, 'raw_cca') & ...
        ~strcmp(allCompsSac, 'baselineVsComp'));
end

% Identify all unique comparisons for attention task
if ~isempty(attnEvents)
    firstAttnEvent = attnEvents{1};
    allCompsAttn = fieldnames(aggSC.(firstAttnEvent));
    allComps.attn = allCompsAttn(startsWith(allCompsAttn, 'raw_') & ...
        ~strcmp(allCompsAttn, 'raw_cca') & ...
        ~strcmp(allCompsAttn, 'baselineVsComp'));
end

% Determine layout dimensions
nRowsSac = length(allComps.sac);
nRowsAttn = length(allComps.attn);
nRows = max(nRowsSac, nRowsAttn);
nCols = nSacEvents + nAttnEvents;

% Get screen dimensions for figure
screenDims = get(0, 'ScreenSize');
screenDims = screenDims(3:4);
figDims = [0.9, 0.5] .* screenDims;
figPos = [screenDims(1:2)/2 - figDims/2, figDims];

% Create figure
fig = figure('Color', 'w', 'Position', figPos, ...
    'MenuBar', 'None', 'ToolBar', 'None');

% Define parameters for mySubPlot
totalWidth = 0.475;        % Total width occupied by all subplots
totalHeight = 0.815;       % Total height occupied by all subplots
leftMargin = 0.025;        % Space on the left side
bottomMargin = 0.065;      % Space at the bottom
horizontalGap = 0.005;     % Horizontal gap between subplots
verticalGap = 0.035;      % Vertical gap between subplots

% Create array to store axes handles
allAx = gobjects(nRows, nCols);

% Loop through saccade events and comparisons
for iCol = 1:nSacEvents
    eventName = sacEvents{iCol};
    eventData = aggSC.(eventName);
    
    for iRow = 1:nRowsSac
        compName = allComps.sac{iRow};
        
        % Create subplot using mySubPlot
        ax = mySubPlot([nRows, nSacEvents, (iRow-1)*nSacEvents + iCol], ...
            totalWidth, totalHeight, leftMargin, bottomMargin, ...
            horizontalGap, verticalGap);
        
        allAx(iRow, iCol) = ax;
        
        % Check if data exists for this comparison
        if isfield(eventData, compName) && isfield(eventData.(compName), 'sig')
            sigMatrix = eventData.(compName).sig;
            sigMatrix = -1 * sigMatrix; 
            nNeurons = size(sigMatrix, 1);
            
            % Get time vector
            if isfield(aggSC, 'binCenters') && isfield(aggSC.binCenters, eventName)
                timeVec = aggSC.binCenters.(eventName);
            else
                timeVec = linspace(-0.5, 1.0, size(sigMatrix, 2));
            end
            
            % Compute preference percentages
            rightPref = sum(sigMatrix == 1, 1) / nNeurons * 100;
            leftPref = -sum(sigMatrix == -1, 1) / nNeurons * 100;
            
            % Plot
            hold(ax, 'on');
            h1 = barStairs(timeVec, rightPref);
            if ~isempty(h1) && length(h1) > 1
                delete(h1(2));
            end
            if ~isempty(h1)
                set(h1(1), 'FaceColor', palette(1,:), 'EdgeColor', 'none');
            end
            
            h2 = barStairs(timeVec, leftPref);
            if ~isempty(h2) && length(h2) > 1
                delete(h2(2));
            end
            if ~isempty(h2)
                set(h2(1), 'FaceColor', palette(2,:), 'EdgeColor', 'none');
            end
            
            % Add event marker
            xline(ax, 0, 'k--');
            
            hold(ax, 'off');
            
            % Add title for first row only
            if iRow == 1
                title(ax, strrep(eventName, '_', ' '), 'Interpreter', 'none');
            end
            
            % Add y-label for first column only
            if iCol == 1
                ylabel(ax, getLabelForComparison(compName));
            end
            
            % Add x-label for bottom row only
            if iRow == nRowsSac
                xlabel(ax, 'Time (s)');
            end
        else
            % No data
            text(0.5, 0.5, 'No data', 'Parent', ax, 'HorizontalAlignment', 'center');
        end
    end
    
    % Fill in any remaining tiles in this column
    for iEmptyRow = nRowsSac+1:nRows
        ax = mySubPlot([nRows, nCols, (iEmptyRow-1)*nCols + iCol], ...
            totalWidth, totalHeight, leftMargin, bottomMargin, ...
            horizontalGap, verticalGap);
        allAx(iEmptyRow, iCol) = ax;
        axis(ax, 'off');
    end
end

% Loop through attention events and comparisons
for iCol = 1:nAttnEvents
    eventName = attnEvents{iCol};
    eventData = aggSC.(eventName);
    colIdx = nSacEvents + iCol;
    
    for iRow = 1:nRowsAttn
        compName = allComps.attn{iRow};
        
        % Create subplot using mySubPlot
        ax = mySubPlot([nRows, nAttnEvents, (iRow-1)*nAttnEvents + iCol], ...
            totalWidth, totalHeight, leftMargin+0.5, bottomMargin, ...
            horizontalGap, verticalGap);
        
        allAx(iRow, colIdx) = ax;
        
        % Check if data exists for this comparison
        if isfield(eventData, compName) && isfield(eventData.(compName), 'sig')
            sigMatrix = aggSC.(eventName).(compName).sig;
            sigMatrix = -1 * sigMatrix;  % Invert to match convention
            nNeurons = size(sigMatrix, 1);
            
            % Get time vector
            if isfield(aggSC, 'binCenters') && isfield(aggSC.binCenters, eventName)
                timeVec = aggSC.binCenters.(eventName);
            else
                timeVec = linspace(-0.5, 1.0, size(sigMatrix, 2));
            end
            
            % Compute preference percentages
            rightPref = sum(sigMatrix == 1, 1) / nNeurons * 100;
            leftPref = -sum(sigMatrix == -1, 1) / nNeurons * 100;
            
            % Plot
            hold(ax, 'on');
            h1 = barStairs(timeVec, rightPref);
            if ~isempty(h1) && length(h1) > 1
                delete(h1(2));
            end
            if ~isempty(h1)
                set(h1(1), 'FaceColor', palette(1,:), 'EdgeColor', 'none');
            end
            
            h2 = barStairs(timeVec, leftPref);
            if ~isempty(h2) && length(h2) > 1
                delete(h2(2));
            end
            if ~isempty(h2)
                set(h2(1), 'FaceColor', palette(2,:), 'EdgeColor', 'none');
            end
            
            % Add event marker
            xline(ax, 0, 'k--');
            
            hold(ax, 'off');
            
            % Add title for first row only
            if iRow == 1
                title(ax, strrep(eventName, '_', ' '), 'Interpreter', 'none');
            end
            
            % Add y-label for first attention column only
            if iCol == 1
                ylabel(ax, getLabelForComparison(compName));
            end
            
            % Add x-label for bottom row only
            if iRow == nRowsAttn
                xlabel(ax, 'Time (s)');
            end
        else
            % No data
            text(0.5, 0.5, 'No data', 'Parent', ax, 'HorizontalAlignment', 'center');
        end
    end
    
    % Fill in any remaining tiles in this column
    for iEmptyRow = nRowsAttn+1:nRows
        ax = mySubPlot([nRows, nCols, (iEmptyRow-1)*nCols + colIdx], ...
            totalWidth, totalHeight, leftMargin, bottomMargin, ...
            horizontalGap, verticalGap);
        allAx(iEmptyRow, colIdx) = ax;
        axis(ax, 'off');
    end
end

% Add overall title
sgtitle('SC Population Preference', 'FontWeight', 'Bold');

% Define which columns should have y-tick labels
yLabelCols = 1;
if ~isempty(attnEvents)
    yLabelCols = [yLabelCols, nSacEvents + 1];
end

% Get all valid axes
validAxes = allAx(isgraphics(allAx));

% Remove y-tick labels except for designated columns
for iCol = 1:nCols
    if ~ismember(iCol, yLabelCols)
        colAxes = allAx(:, iCol);
        colAxes = colAxes(isgraphics(colAxes));
        set(colAxes, 'YTickLabel', []);
        
        % Also remove y-labels
        for iAx = 1:length(colAxes)
            set(get(colAxes(iAx), 'YLabel'), 'String', '');
        end
    end
end

% Remove x-tick labels except for bottom row of plots:
set(allAx(1:(nRowsSac-1),1:nSacEvents), 'XTickLabel', []);
set(allAx(1:(nRowsAttn-1),1:nAttnEvents), 'XTickLabel', []);

% Set y-limits based on task type
yQuantPref = 25;

% Handle saccade task plots
for iRow = 1:nRowsSac
    sacAxes = allAx(iRow, 1:nSacEvents);
    sacAxes = sacAxes(isgraphics(sacAxes));
    if ~isempty(sacAxes)
        [~, yLims] = outerLims(sacAxes);
        yLims = quantizeLimits(yLims, yQuantPref);
        set(sacAxes, 'YLim', yLims, 'YTick', ...
            yLims(1):yQuantPref:yLims(2));
    end
end

% Handle attention task plots
for iRow = 1:nRowsAttn
    attnCols = (nSacEvents+1):(nSacEvents+nAttnEvents);
    attnAxes = allAx(iRow, attnCols);
    attnAxes = attnAxes(isgraphics(attnAxes));
    if ~isempty(attnAxes)
        [~, yLims] = outerLims(attnAxes);
        yLims = quantizeLimits(yLims, yQuantPref);
        set(attnAxes, 'YLim', yLims, 'YTick', ...
            yLims(1):yQuantPref:yLims(2));
    end
end

% Add legend to top-right plot
if isgraphics(allAx(1,1))
    legend(allAx(1,1), {'Ipsi', 'Contra'}, 'Location', 'best', ...
        'Box','Off');
end

% Set common formatting
set(validAxes, 'TickDir', 'Out', 'LineWidth', 0.75, 'XGrid', 'Off', ...
    'YGrid', 'Off', 'Box', 'off');

% Set common x-limits within columns
xQuant = 0.1;
for iCol = 1:nCols
    colAxes = allAx(:, iCol);
    colAxes = colAxes(isgraphics(colAxes));
    if ~isempty(colAxes)
        xLims = quantizeLimits(outerLims(colAxes), xQuant);
        set(colAxes, 'XLim', xLims);
    end
end

% Save the figure
if nargin >= 2
    figName = fullfile(saveDir, 'sc_population_preference');
    pdfSave(figName, fig.Position(3:4)/72, fig);
    giveFeed('Saved figure');
end
end

function label = getLabelForComparison(compName)
% Returns a readable label for a comparison name
compName = strrep(compName, 'raw_', '');

switch compName
    case 'targetLoc'
        label = 'Target Location';
    case 'saliencyTgtContra'
        label = 'Target Saliency';
    case 'rewardTgtContra'
        label = 'Target Reward';
    case 'cueLoc'
        label = 'Cue Location';
    case 'cueChgLoc'
        label = 'Cue Change Location';
    case 'cueChgInHitVsMiss'
        label = 'Hit vs Miss';
    case 'saliency'
        label = 'Saliency';
    case 'reward'
        label = 'Reward';
    otherwise
        label = strrep(compName, '_', ' ');
end
end

function newLims = quantizeLimits(limits, q)
% Quantize axis limits to the nearest multiple of q
lowerLimit = limits(1);
upperLimit = limits(2);
newLowerLimit = floor(lowerLimit / q) * q;
newUpperLimit = ceil(upperLimit / q) * q;
newLims = [newLowerLimit, newUpperLimit];
end