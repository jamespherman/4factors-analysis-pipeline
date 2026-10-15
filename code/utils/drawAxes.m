function h = drawAxes(ax, X, Y, varargin)
%DRAWAXES Create disconnected L-shaped axes using annotation objects.
%
%   h = drawAxes(ax, X, Y)
%   h = drawAxes(ax, X, Y, tickLengthProp)
%   h = drawAxes(ax, X, Y, tickLengthProp, nXTicks)
%   h = drawAxes(ax, X, Y, tickLengthProp, nXTicks, nYTicks)
%
%   MATLAB's built-in axes always form a closed box or connected L-shape.
%   This function creates open, disconnected axes where the X and Y axis
%   lines do not touch at the origin, which is useful for publication-
%   quality figures (e.g., ROC plots, scatter plots).
%
%   INPUTS
%       ax              - Handle to target axes object
%       X               - 3-element vector [xAxisPos, xMin, xMax] in data units
%                           xAxisPos: x-coordinate where Y-axis line is drawn
%                           xMin:     x-coordinate where X-axis line begins
%                           xMax:     x-coordinate where X-axis line ends
%       Y               - 3-element vector [yAxisPos, yMin, yMax] in data units
%                           yAxisPos: y-coordinate where X-axis line is drawn
%                           yMin:     y-coordinate where Y-axis line begins
%                           yMax:     y-coordinate where Y-axis line ends
%       tickLengthProp  - (Optional) Tick length as proportion of axes 
%                         height/width. Default: 0.035
%       nXTicks         - (Optional) Number of ticks on X-axis. Default: 5
%       nYTicks         - (Optional) Number of ticks on Y-axis. Default: 5
%
%   OUTPUT
%       h               - Column vector of annotation line handles for all
%                         created axis lines and ticks (tagged 'X Axis', 
%                         'Y Axis', 'X Tick', 'Y Tick')
%
%   CREATING DISCONNECTED AXES
%       To create a gap at the origin, set axis positions outside the data
%       range:
%           xAxisPos < xMin  →  Y-axis is left of where X-axis starts
%           yAxisPos < yMin  →  X-axis is below where Y-axis starts
%
%   EXAMPLES
%       % Connected axes (standard L-shape, no gap):
%       X = [0, 0, 1];    Y = [0, 0, 1];
%       drawAxes(ax, X, Y);
%
%       % Disconnected axes (gap at origin):
%       X = [-0.08, 0, 1];    Y = [-0.08, 0, 1];
%       drawAxes(ax, X, Y);
%
%       % Custom ticks (longer ticks, 3 on X, 7 on Y):
%       drawAxes(ax, X, Y, 0.05, 3, 7);
%
%       % Draw only X-axis (set ax.YColor='none' before calling):
%       ax.YColor = 'none';
%       drawAxes(ax, X, Y);
%
%   NOTES
%       - Sets ax.XLim to [xMin, xMax] and ax.YLim to [yMin, yMax]
%       - Hides MATLAB's built-in axes after drawing (sets XColor/YColor 
%         to 'none')
%       - Axes are drawn as annotation objects on the parent figure
%       - To skip drawing an axis, set its color to 'none' before calling
%
%   See also ANNOTATION, AXES

% -------------------------------------------------------------------------
% Input parsing
% -------------------------------------------------------------------------

% Support legacy 4-element X/Y vectors (4th element is unused)
if numel(X) > 3, X = X(1:3); end
if numel(Y) > 3, Y = Y(1:3); end

% Parse optional arguments
tickLengthProp = 0.035;
nXTicks = 5;
nYTicks = 5;

if nargin >= 4, tickLengthProp = varargin{1}; end
if nargin >= 5, nXTicks = varargin{2}; end
if nargin >= 6, nYTicks = varargin{3}; end

% -------------------------------------------------------------------------
% Setup
% -------------------------------------------------------------------------

figParent = get(ax, 'Parent');

% Store and normalize units for coordinate calculations
prevUnits = get(ax, 'Units');
set(ax, 'Units', 'Normalized');

% Check which axes should be drawn (skip if color is already 'none')
drawXAxis = ~ischar(ax.XColor);
drawYAxis = ~ischar(ax.YColor);

% Set axis limits from input vectors
ax.XLim = X([2 3]);
ax.YLim = Y([2 3]);
drawnow;

% -------------------------------------------------------------------------
% Coordinate conversion
% -------------------------------------------------------------------------

% Compute tick locations in data units
xTickData = linspace(X(2), X(3), nXTicks);
yTickData = linspace(Y(2), Y(3), nYTicks);

% Convert all coordinates from data units to normalized figure units
[nX, nY] = data2norm(ax, X, Y);
[xTickNorm, yTickNorm] = data2norm(ax, xTickData(:), yTickData(:));

% Calculate tick lengths in normalized units (proportional to axes size)
axPos = get(ax, 'Position');
xTickLen = axPos(4) * tickLengthProp;  % X-ticks: proportional to height
yTickLen = axPos(3) * tickLengthProp;  % Y-ticks: proportional to width

% -------------------------------------------------------------------------
% Draw axes and ticks
% -------------------------------------------------------------------------

h = [];  % Accumulate handles for output

% Draw X-axis line and ticks
if drawXAxis
    % X-axis line: horizontal line at yAxisPos, from xMin to xMax
    hXLine = annotation(figParent, 'line', nX(2:3), nY(1)*[1 1], 'Color', 'k');
    hXLine.Tag = 'X Axis';
    h = [h; hXLine];
    
    % X-axis ticks: vertical lines extending downward from axis
    if nXTicks > 0
        xTickYCoords = [nY(1), nY(1) - xTickLen];
        for i = 1:nXTicks
            hTick = annotation(figParent, 'line', xTickNorm(i)*[1 1], xTickYCoords, 'Color', 'k');
            hTick.Tag = 'X Tick';
            h = [h; hTick];
        end
    end
end

% Draw Y-axis line and ticks
if drawYAxis
    % Y-axis line: vertical line at xAxisPos, from yMin to yMax
    hYLine = annotation(figParent, 'line', nX(1)*[1 1], nY(2:3), 'Color', 'k');
    hYLine.Tag = 'Y Axis';
    h = [h; hYLine];
    
    % Y-axis ticks: horizontal lines extending leftward from axis
    if nYTicks > 0
        yTickXCoords = [nX(1), nX(1) - yTickLen];
        for i = 1:nYTicks
            hTick = annotation(figParent, 'line', yTickXCoords, yTickNorm(i)*[1 1], 'Color', 'k');
            hTick.Tag = 'Y Tick';
            h = [h; hTick];
        end
    end
end

% -------------------------------------------------------------------------
% Cleanup
% -------------------------------------------------------------------------

% Hide MATLAB's built-in axes visuals
ax.PositionConstraint = 'innerposition';
ax.XColor = 'none';
ax.YColor = 'none';

% Restore original units
set(ax, 'Units', prevUnits);

end


function [normX, normY] = data2norm(ax, dataX, dataY)
%DATA2NORM Convert data coordinates to normalized figure coordinates.
%
%   [normX, normY] = data2norm(ax, dataX, dataY)
%
%   Applies a linear transformation from the axes' data coordinate system
%   to the parent figure's normalized [0,1] coordinate system.

% Linear transformation: norm = m * data + b
mX = ax.Position(3) / diff(ax.XLim);
mY = ax.Position(4) / diff(ax.YLim);
bX = ax.Position(1) - mX * ax.XLim(1);
bY = ax.Position(2) - mY * ax.YLim(1);

normX = mX * dataX + bX;
normY = mY * dataY + bY;

end