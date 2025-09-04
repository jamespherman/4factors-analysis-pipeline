function varargout = alignAndBinSpikes(varargin)
%
% [alignedSpikes, binnedCountsArray, binCenters] = alignAndBinSpikes(...
%   spikeTimes, eventTimes, minTime, maxTime, binWidth, [slidingBinWidth],
%   [slideTime]).
%
% alignAndBinSpikes:    Aligns spike times to given event times and bins 
%                       them to create a histogram.

% INPUTS:
%   spikeTimes: Vector of spike times for a single spikeTimes.
%   eventTimes: Times at which the events occur.
%   minTime, maxTime: Define the time window for alignment.
%   binWidth: Width of the bins for the histogram.
%   nEventTimes: Number of event times.
%   nBins: Number of bins.
%   slidingBinWidth: Width of SLIDING bins (seconds).
%   slideTime: How much to slide each bin by (seconds).

% OUTPUTS:
%   alignedSpikes:      Cell array with spikes aligned to each event time.
%   binnedCountsArray:  Matrix with the histogram counts for spikes aligned 
%                       to each event time.
%   eventNumberVector:  A vector of integers the same length as
%                       "alignedSpikes" indicating which event each spike 
%                       was aligned to.

% use input parser:
p = inputParser;
addRequired(p, 'spikeTimes');
addRequired(p, 'eventTimes');
addRequired(p, 'minTime');
addRequired(p, 'maxTime');
addRequired(p, 'binWidth');
addOptional(p, 'slidingBinWidth', 0.2);
addOptional(p, 'slideTime', 0.001);

parse(p, varargin{:});
spikeTimes      = p.Results.spikeTimes;
eventTimes      = p.Results.eventTimes;
minTime         = p.Results.minTime;
maxTime         = p.Results.maxTime;
binWidth        = p.Results.binWidth;
slidingBinWidth = p.Results.slidingBinWidth;
slideTime       = p.Results.slideTime;

% define number of output arguments and argument names
nOutputs = nargout;
varargout = cell(1,nOutputs);
argOutNames = {'alignedSpikes', 'binnedCountsArray', ...
    'binCenters', 'eventNumberVector', 'slidingBinnedCountsArray', ...
    'slidingBinTimes'};

% Define values and make variables to store outputs
binStarts = (round(minTime/0.001)*0.001):binWidth:(round((maxTime - ...
    binWidth)/0.001)*0.001);
binEnds   = binStarts + binWidth;
binCenters = binStarts + (binWidth/2);
nBins = length(binStarts);
nEventTimes = length(eventTimes);
nSlidingBins = round((maxTime - minTime)*1000) - ...
    round(slideTime * 1000);
alignedSpikes = cell(nEventTimes, 1);
binnedCountsArray = zeros(nEventTimes, nBins);
eventNumberVector = alignedSpikes;
slidingBinTimes = zeros(nSlidingBins, 1);
slidingBinnedCountsArray = zeros(nEventTimes, nSlidingBins);

% loop over event times
for i = 1:nEventTimes

    % logically index spike times that fall within our window of interest
    % relative to the current event time:
    g = spikeTimes > eventTimes(i) + minTime & spikeTimes < ...
        eventTimes(i) + maxTime;

    % logically index spike times that fall within our window of interest
    % +/- an additional half a slidingBinWidth so we end up with only full
    % sized sliding bin widths of interest:
    gSl = spikeTimes > eventTimes(i) + minTime - slidingBinWidth/2 ...
        & spikeTimes < eventTimes(i) + maxTime + slidingBinWidth/2;

    % Retreive and align spikes in our window:
    tempAlignedSpikes = spikeTimes(g) - eventTimes(i);

    % store aligned spikes in cell array
    alignedSpikes{i} = tempAlignedSpikes;

    % store binned counts
    %     binnedCountsArray(i,:) = histcounts(tempAlignedSpikes, ...
    %         'BinLimits',[minTime, maxTime],'BinWidth',binWidth);
    for j = 1:nBins
        binnedCountsArray(i, j) = nnz(tempAlignedSpikes > binStarts(j) ...
            & tempAlignedSpikes < binEnds(j));
    end

    % store sliding binned counts (if desired):
    if nargout > 3
        try
        [counts, sBinTimes] = slidingBinnedCounts(...
            slidingBinWidth, spikeTimes(gSl) - eventTimes(i), slideTime,...
            minTime - slidingBinWidth/2, minTime - slidingBinWidth/2, ...
            maxTime + slidingBinWidth/2);
        g2 = round(sBinTimes, 3) > round(minTime, 3) & ...
            round(sBinTimes, 3) < round(maxTime, 3);
        slidingBinnedCountsArray(i, :) = counts(g2);
        if i == 1
            slidingBinTimes = sBinTimes(g2);
        end
        catch me
            keyboard
        end
    end

    % Store "trial number" (event number):
    eventNumberVector{i} = i * ones(nnz(g), 1);
end

% logically index the sliding bin centers that fall outside of our
% window of interest:
gOut = slidingBinTimes < minTime | slidingBinTimes > maxTime;

% get rid of slidingBinTimes outside our window of interest, and
% slidingBinnedCountsArray columns outside our window of interest.
if nargout > 3
    slidingBinTimes(gOut) = [];
    slidingBinnedCountsArray(:, gOut) = [];
end

% concatenate
alignedSpikes = vertcat(alignedSpikes{:});
eventNumberVector = vertcat(eventNumberVector{:});

% return outputs
for i = 1:nargout
    varargout{i} = eval(argOutNames{i});
end
