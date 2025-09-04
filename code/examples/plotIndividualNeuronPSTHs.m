function plotIndividualNeuronPSTHs(data, sessionName, saveFolder)
% plotIndividualNeuronPSTHs generates a tiled figure for each neuron,
% plotting only PSTHs for each event and each comparison.
%
%   plotIndividualNeuronPSTHs(data, sessionName, saveFolder)
%
%   data        - structure containing analysisConfig and rawArray data
%   sessionName - session file name (will be cleaned for display and file naming)
%   saveFolder  - folder where the PDF files will be saved

useRaw = true;

% Clean up sessionName for display and file naming
sessionNameClean = strrep(strrep(sessionName, '_', ' '), '.mat', '');
sessionNameFile  = strrep(sessionNameClean, ' ', '_');

% Get all event names from analysisConfig
events = fieldnames(data.analysisConfig);

% Separate events by task type
sacEvents = events(startsWith(events, 'sac_'));
attnEvents = events(startsWith(events, 'attn_'));
nSacEvents = length(sacEvents);
nAttnEvents = length(attnEvents);

% Combine events in order (first saccade events, then attention events)
allEvents = [sacEvents; attnEvents];
nEvents = length(allEvents);

% Determine which columns should have y-axis labels
yLabelCols = [1, nSacEvents + 1]; % First column of each task type

% Determine maximum number of comparisons across all events
maxComps = 0;
for iEvent = 1:nEvents
    cfg = data.analysisConfig.(allEvents{iEvent});
    [comps, ~] = findComps(cfg);
    if useRaw
        comps = comps(startsWith(comps, 'raw_'));
    else
        comps = comps(~startsWith(comps, 'raw_'));
    end
    maxComps = max(maxComps, length(comps));
end

% Determine the maximum number of neurons across all events
nNeurons = 0;
for iEvent = 1:nEvents
    currNeurons = size(data.analysisConfig.(allEvents{iEvent}).rawArray, 2);
    nNeurons = max(nNeurons, currNeurons);
end

for neuronIdx = 1:nNeurons
    % Set up figure and tiled layout (rows = maxComps, columns = nEvents)
    screenDims = get(0, 'ScreenSize');
    screenDims = screenDims(3:4);
    figDims = [0.6 0.6] .* screenDims;
    figPos  = (screenDims - figDims) / 2;
    figPos  = [figPos, figDims];
    
    fig = figure('Color', 'w', 'Position', figPos, ...
        'MenuBar', 'None', 'ToolBar', 'None');
    tl = tiledlayout(maxComps, nEvents, 'TileSpacing', 'Compact', ...
        'Padding', 'Compact');
    
    % Load color palette (using richColors function)
    palette = richColors;
    palette = palette([6, 12], :);
    
    % Preallocate an axes handle array (for later adjustment of limits)
    axHandles = gobjects(maxComps, nEvents);
    
    % Loop over each event (each column)
    for iEvent = 1:nEvents
        eventName = allEvents{iEvent};
        cfg = data.analysisConfig.(eventName);
        
        % Get comparison fields for this event
        [comps, ~] = findComps(cfg);
        if useRaw
            comps = comps(startsWith(comps, 'raw_'));
        else
            comps = comps(~startsWith(comps, 'raw_'));
        end
        nComps = length(comps);
        
        % Determine if this is a saccade or attention task and get appropriate comparisons
        if contains(eventName, 'sac_')
            compToPass = data.comparisons.gSac;
            taskType = 'saccade';
        else
            compToPass = data.comparisons.attn;
            taskType = 'attention';
        end
        
        % Loop over each comparison (each row)
        for iComp = 1:nComps
            compName = comps{iComp};
            
            % Get condition vectors and labels for this comparison using the new approach
            [cond1, cond2, label1, label2] = getConditionInfo(cfg, compName, compToPass);
            
            % Compute tile index (tiles are filled row-by-row)
            tileIdx = (iComp - 1) * nEvents + iEvent;
            ax = nexttile(tl, tileIdx);
            axHandles(iComp, iEvent) = ax;
            
            % Check if this event's rawArray has the current neuron index
            if neuronIdx <= size(cfg.rawArray, 2)
                % Compute the PSTH for this neuron and conditions
                psth1 = squeeze(mean(cfg.rawArray(cond1, neuronIdx, :), 1));
                psth2 = squeeze(mean(cfg.rawArray(cond2, neuronIdx, :), 1));
                timeVec = cfg.binCenters;
                
                % Plot the PSTHs using barStairs (custom plotting function)
                h1 = barStairs(timeVec, psth1, false);
                set(h1, 'Color', palette(1,:));
                hold on;
                h2 = barStairs(timeVec, psth2, false);
                set(h2, 'Color', palette(2,:));
                
                % Add event marker
                xline(0, 'k--');
            else
                % If the neuron index exceeds the number available for this event,
                % leave the tile empty and add a note.
                axis off;
                text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center', 'FontSize', 10);
            end
            
            % Add title only for the first row and remove colon and subsequent text
            if iComp == 1
                colonIdx = strfind(eventName, ':');
                if ~isempty(colonIdx)
                    titleStr = eventName(1:colonIdx(1)-1);
                else
                    titleStr = strrep(eventName, '_', ' ');
                end
                title(titleStr, 'FontSize', 10);
            end
            
            % Only include x-axis tick labels for the final row
            if iComp == maxComps || (iComp == 3 && iEvent > 5)
                xlabel('Time (s)');
            else
                set(gca, 'XTickLabel', []);
            end
            
            % Only include y-axis tick labels for the first column of each task type
            if ~ismember(iEvent, yLabelCols)
                set(gca, 'YTickLabel', []);
            else
                ylabel('Mean Spike Count');
            end
            
            % Add legends 
            % For first column of each task type, add legend for each comparison
            if ismember(iEvent, yLabelCols)
                if ishandle(h1) && ishandle(h2)
                    % Customize labels based on task type and comparison
                    if strcmp(taskType, 'saccade')
                        switch compName
                            case 'raw_targetLoc' 
                                leg1 = 'Target Contra';
                                leg2 = 'Target Ipsi';
                            case {'raw_saliency', 'raw_saliencyTgtContra'}
                                leg1 = 'High Salience';
                                leg2 = 'Low Salience';
                            case {'raw_reward', 'raw_rewardTgtContra'}
                                leg1 = 'High Reward';
                                leg2 = 'Low Reward';
                            otherwise
                                keyboard
                                leg1 = label1;
                                leg2 = label2;
                        end
                    else % attention task
                        switch compName
                            case {'raw_cueLoc', 'raw_cueLocation'}
                                leg1 = 'Cue in RF';
                                leg2 = 'Cue out RF';
                            case {'raw_cueChgInHitVsMiss', ...
                                    'raw_hitMissCueIn'}
                                leg1 = 'Hit in RF';
                                leg2 = 'Miss in RF';
                            case 'raw_cueChgLoc'
                                leg1 = 'Cue Change in RF';
                                leg2 = 'Cue Change out RF';
                            otherwise
                                keyboard
                                leg1 = label1;
                                leg2 = label2;
                        end
                    end
                    legend([h1, h2], {leg1, leg2}, 'Location', 'northeast');
                end
            end
        end
        
        % Fill empty tiles if an event has fewer comparisons than maxComps
        for iComp = nComps+1:maxComps
            tileIdx = (iComp - 1) * nEvents + iEvent;
            nexttile(tl, tileIdx);
            axis off;
        end
    end
    
    % --- Uniform y-axis limits adjustment ---
    yQuantPref = 0.25;
    
    % Saccade task plots (columns 1:nSacEvents)
    sacCols = 1:nSacEvents;
    sacAxes = axHandles(:, sacCols);
    sacAxes = sacAxes(arrayfun(@(ax) isgraphics(ax), sacAxes));
    if ~isempty(sacAxes)
        yLimsSac = outerLims(sacAxes);
        yLimsSac = quantizeLimits(yLimsSac, yQuantPref);
        for ax = sacAxes'
            set(ax, 'YLim', yLimsSac, 'YTick', yLimsSac(1):yQuantPref:yLimsSac(2));
        end
    end

    % Attention task plots (columns nSacEvents+1:nEvents)
    attnCols = (nSacEvents+1):nEvents;
    attnAxes = axHandles(:, attnCols);
    attnAxes = attnAxes(arrayfun(@(ax) isgraphics(ax), attnAxes));
    if ~isempty(attnAxes)
        yLimsAttn = outerLims(attnAxes);
        yLimsAttn = quantizeLimits(yLimsAttn, yQuantPref);
        for ax = attnAxes'
            set(ax, 'YLim', yLimsAttn, 'YTick', yLimsAttn(1):yQuantPref:yLimsAttn(2));
        end
    end
    % --- End y-axis adjustment ---
    
    % Add an overall title for the figure
    sgtitle(tl, sprintf('%s - Neuron %d', sessionNameClean, neuronIdx));
    set([sacAxes; attnAxes], 'TickDir', 'Out', 'LineWidth', 0.75, ...
        'XColor', 'k', 'YColor', 'k')
    delete(findall(0, 'Type', 'Legend'))

    % Save the figure as a PDF file
    fileName = sprintf('%s_neuron%d.pdf', sessionNameFile, neuronIdx);
    % Use standard print function if savem is not available
    if exist('savem', 'file')
        savem(saveFolder, fileName, fig);
    else
        fullPath = fullfile(saveFolder, fileName);
        print(fig, fullPath, '-dpdf', '-r300');
    end
    close(fig);
end
end

function [cond1, cond2, label1, label2] = getConditionInfo(cfg, compName, compToPass)
% [cond1, cond2, label1, label2] = getConditionInfo(cfg, compName, compToPass)
%
%   Retrieves condition vectors and labels for a given comparison
%   name from the configuration structure.
%
%   INPUTS:
%     cfg             : Configuration structure containing condition
%                       variables in cfg.condVars.
%     compName        : Name of the comparison (string). Can include
%                       'raw_' prefix, which will be removed.
%     compToPass      : Structure array defining comparison names,
%                       conditions, and labels.
%
%   OUTPUTS:
%     cond1           : Logical vector for condition 1.
%     cond2           : Logical vector for condition 2.
%     label1          : Label for condition 1 (string).
%     label2          : Label for condition 2 (string).

% Remove 'raw_' prefix if it exists.
if startsWith(compName, 'raw_')
    compName = extractAfter(compName, 4);
end

% Look up the comparison info from the provided comparisons structure.
idx = find(strcmpi({compToPass.name}, compName), 1);
if isempty(idx)
    error(['Comparison information not found for: ', compName]);
end
cmp = compToPass(idx);

% Extract condition vectors from cfg.condVars using the names from cmp.
cond1 = cfg.condVars.(cmp.condition1);
cond2 = cfg.condVars.(cmp.condition2);

% Use the condition names as labels (remove 'is' and 'Trial' to make them more readable)
label1 = regexprep(cmp.condition1, 'is|Trial', '');
label2 = regexprep(cmp.condition2, 'is|Trial', '');
end

function [isComp, notComp] = findComps(cfg)
% Get all field names
fNames = fieldnames(cfg);
% Determine which fields are structures
isStruct = structfun(@isstruct, cfg);
% Check for 'singleNeuron' subfield in each structure field
hasSingleNeuron = false(size(isStruct));
for i = 1:length(fNames)
    if isStruct(i)
        hasSingleNeuron(i) = isfield(cfg.(fNames{i}), 'singleNeuron');
    end
end
% Only consider fields that are structures
validFields = isStruct;
% Get field names that are structures with and without the subfield
isComp = fNames(validFields & hasSingleNeuron);
notComp = fNames(validFields & ~hasSingleNeuron);
end

function yLims = outerLims(axArray)
% Compute combined y-limits for an array of axes.
ymin = inf;
ymax = -inf;
for k = 1:length(axArray)
    if isgraphics(axArray(k))
        yl = get(axArray(k), 'YLim');
        ymin = min(ymin, yl(1));
        ymax = max(ymax, yl(2));
    end
end
yLims = [ymin, ymax];
end

function newLims = quantizeLimits(limits, q)
lowerLimit = limits(1);
upperLimit = limits(2);
newLowerLimit = floor(lowerLimit / q) * q;
newUpperLimit = ceil(upperLimit / q) * q;
newLims = [newLowerLimit, newUpperLimit];
end