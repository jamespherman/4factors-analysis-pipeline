function h = mySubPlot(A, varargin)
% MYSUBPLOT Create custom subplot with precise positioning control
%
%   h = mySubPlot(A)
%   h = mySubPlot(A, 'PropertyName', PropertyValue, ...)
%
% -------------------------------------------------------------------------
% CORE CONCEPT
% -------------------------------------------------------------------------
%   Each call to mySubPlot is INDEPENDENT. The function calculates a
%   position based on: "If this figure were divided into an [R x C] grid,
%   where would slot #idx be?" Different calls can use different grids—
%   they do not interact or check for overlap. This enables complex layouts
%   by mixing grid definitions across calls.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   A : A 1x3 vector specifying [Rows, Columns, Index]
%       - Rows and Columns can be fractional (see FRACTIONAL GRIDS below)
%       - Index fills row-wise, starting at top-left
%
% OUTPUT:
%   h : Handle to the created axes. The axes becomes current, so subsequent
%       plot commands will draw into it:
%           ax1 = mySubPlot([2,1,1]);
%           plot(x, y);               % Draws into ax1
%           ax2 = mySubPlot([2,1,2]);
%           plot(ax2, x, z);          % Explicit handle also works
%
% -------------------------------------------------------------------------
% INDEX ORDERING (row-wise, top-to-bottom)
% -------------------------------------------------------------------------
%   For A = [2, 3, idx]:
%
%       +---+---+---+
%       | 1 | 2 | 3 |   <- Row 1 (TOP)
%       +---+---+---+
%       | 4 | 5 | 6 |   <- Row 2 (BOTTOM)
%       +---+---+---+
%
%   For A = [3, 2, idx]:
%
%       +---+---+
%       | 1 | 2 |   <- Row 1 (TOP)
%       +---+---+
%       | 3 | 4 |   <- Row 2
%       +---+---+
%       | 5 | 6 |   <- Row 3 (BOTTOM)
%       +---+---+
%
% -------------------------------------------------------------------------
% OPTIONAL PARAMETERS (Name-Value Pairs)
% -------------------------------------------------------------------------
%   'Width'        : Total width of subplot area (0-1, default 0.85)
%   'Height'       : Total height of subplot area (0-1, default 0.85)
%   'LeftMargin'   : Left offset from figure edge (0-1, default 0.07)
%   'BottomMargin' : Bottom offset from figure edge (0-1, default 0.07)
%   'WidthSpacing' : Horizontal gap between columns (0-1, default 0.015)
%   'HeightSpacing': Vertical gap between rows (0-1, default 0.035)
%
% -------------------------------------------------------------------------
% QUICK REFERENCE - COMMON LAYOUTS
% -------------------------------------------------------------------------
%   2 rows, equal height:           [2,1,1], [2,1,2]
%   3 rows, equal height:           [3,1,1], [3,1,2], [3,1,3]
%   2x2 grid:                       [2,2,1], [2,2,2], [2,2,3], [2,2,4]
%   Top 2/3, bottom 1/3:            [3/2,1,1], [3,1,3]
%   Top 3/4, bottom 1/4:            [4/3,1,1], [4,1,4]
%   Left 1/3, right 2/3:            [1,3,1],   [1,3/2,3]  % Note: idx 3
%   Left 2/3, right 1/3:            [1,3/2,1], [1,3,3]
%
%   LIMITATION: Fractional grids only work for FIRST panels (top/left).
%   For "Bottom 3/4" or "Right 2/3" as a SINGLE large panel, the ceil()
%   logic prevents proper positioning. Instead, stack standard grid panels:
%       Top 1/4, bottom 3/4:        [4,1,1], then [4,1,2],[4,1,3],[4,1,4]
%       Left 1/4, right 3/4:        [4,1,1], then [1,4,2],[1,4,3],[1,4,4]
%
% -------------------------------------------------------------------------
% PATTERN 1: FRACTIONAL GRIDS (The "Inverse Rule")
% -------------------------------------------------------------------------
%   To create a panel occupying fraction P of a dimension, set the grid
%   value to 1/P.
%
%   FORMULA:  GridValue = 1 / DesiredFraction
%
%   IMPORTANT: Fractional grids only work reliably with INDEX 1 (top-left).
%   For other positions, use integer grids or viewport partitioning.
%
%   Examples:
%     - Want 75% height?  1/0.75 = 4/3   -> use [4/3, C, 1]
%     - Want 50% height?  1/0.50 = 2     -> use [2, C, 1]
%     - Want 33% height?  1/0.33 = 3     -> use [3, C, 1]
%     - Want 40% width?   1/0.40 = 2.5   -> use [R, 2.5, 1]
%
%   COMPLETE EXAMPLE: Top 75%, bottom 25% split into two panels
%
%       figure('Color','w');
%       mySubPlot([4/3, 1, 1]);       % Top 75%: row 1 of a 1.33-row grid
%       title('Large top panel');
%
%       mySubPlot([4, 2, 7]);         % Bottom-left 25%: row 4 of 4, col 1
%       title('Small bottom-left');
%
%       mySubPlot([4, 2, 8]);         % Bottom-right 25%: row 4 of 4, col 2
%       title('Small bottom-right');
%
%   WHY THIS WORKS:
%     - [4/3,1,1] creates a panel 1/(4/3) = 75% tall, positioned at top
%     - [4,2,7] and [4,2,8] use a 4x2 grid; indices 7,8 are the bottom row
%     - The 4-row grid's bottom row is exactly 25% of figure height
%
% -------------------------------------------------------------------------
% PATTERN 2: MIXED GRID RESOLUTIONS
% -------------------------------------------------------------------------
%   Instead of "spanning" multiple cells (which this function does NOT
%   support), use different grid definitions for different panel sizes.
%
%   EXAMPLE: Top half as one panel, bottom half as four panels
%
%       figure('Color','w');
%       mySubPlot([2,1,1]);           % Top half (2-row grid, index 1)
%       title('Wide top');
%
%       mySubPlot([4,2,5]);           % Bottom-left quadrant of bottom half
%       mySubPlot([4,2,6]);           % Bottom-right quadrant of bottom half  
%       mySubPlot([4,2,7]);           % etc.
%       mySubPlot([4,2,8]);
%
%   EXAMPLE: Left column narrow, right column wide
%   For layouts mixing fractional widths with multiple rows, use Pattern 3
%   (viewport partitioning) instead—it avoids the ceil() index limitation.
%
%       figure('Color','w');
%       % Left 1/3: use viewport partitioning
%       opts_left = {'Width', 0.28, 'LeftMargin', 0.07};
%       mySubPlot([2,1,1], opts_left{:});  title('Left top');
%       mySubPlot([2,1,2], opts_left{:});  title('Left bottom');
%
%       % Right 2/3: offset viewport
%       opts_right = {'Width', 0.56, 'LeftMargin', 0.38};
%       mySubPlot([2,1,1], opts_right{:}); title('Right top');
%       mySubPlot([2,1,2], opts_right{:}); title('Right bottom');
%
% -------------------------------------------------------------------------
% PATTERN 3: INDEPENDENT COLUMN GROUPS (Viewport Partitioning)
% -------------------------------------------------------------------------
%   Use 'LeftMargin' and 'Width' to create separate vertical regions with
%   independent grid densities.
%
%   EXAMPLE: Three independent column regions
%
%       figure('Color','w', 'Position', [100 100 1200 600]);
%
%       % Column region 1 (left): 4 stacked panels
%       opts1 = {'Width', 0.25, 'LeftMargin', 0.07};
%       for i = 1:4
%           mySubPlot([4,1,i], opts1{:});
%           title(sprintf('Col1 Panel %d', i));
%       end
%
%       % Column region 2 (middle): 2x2 grid
%       opts2 = {'Width', 0.30, 'LeftMargin', 0.40};
%       for i = 1:4
%           mySubPlot([2,2,i], opts2{:});
%           title(sprintf('Col2 Panel %d', i));
%       end
%
%       % Column region 3 (right): 3 stacked panels
%       opts3 = {'Width', 0.20, 'LeftMargin', 0.75};
%       for i = 1:3
%           mySubPlot([3,1,i], opts3{:});
%           title(sprintf('Col3 Panel %d', i));
%       end
%
% -------------------------------------------------------------------------
% COMMON MISTAKES (GOTCHAS)
% -------------------------------------------------------------------------
%   1. DO NOT use vector indices like subplot does:
%        WRONG:  mySubPlot([2, 3, [1,2,3]])    % Will error
%        RIGHT:  mySubPlot([2, 1, 1])          % Change grid instead
%
%   2. Fractional grids ONLY work for top/left panels (index 1). Due to
%      ceil() in the position math, fractional grids with index > 1 will
%      place panels off-screen. For large bottom/right panels, use
%      multiple standard grid calls instead.
%        WRONG:  mySubPlot([4/3, 1, 2])        % Off-screen (bottom 3/4)
%        RIGHT:  mySubPlot([4,1,2]), mySubPlot([4,1,3]), mySubPlot([4,1,4])
%
%   3. Panels CAN overlap if your grid definitions conflict. This is
%      by design—the function does not track previous calls.
%
%   4. Index 1 is always TOP-LEFT, not bottom-left. This differs from
%      some plotting conventions where y=0 is at bottom.
%
%   5. The function creates a new axes every call. It does not reuse or
%      clear existing axes at that position.
%
% -------------------------------------------------------------------------
% Original by jph - 1/23/2015
% Modernized with inputParser - 2/25/2025
% Documentation updated for LLM context - 1/10/2026
% -------------------------------------------------------------------------

% === IMPLEMENTATION ===

% Create input parser
p = inputParser;
p.FunctionName = 'mySubPlot';
p.CaseSensitive = false;
p.KeepUnmatched = true;

% Validate subplot grid specification
validateA = @(x) validateattributes(x, {'numeric'}, ...
    {'vector', 'numel', 3, 'positive'});
p.addRequired('A', validateA);

% Add optional parameters with validation
p.addParameter('Width', 0.85, @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'positive', '<=', 1}));
p.addParameter('Height', 0.85, @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'positive', '<=', 1}));
p.addParameter('LeftMargin', 0.07, @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'nonnegative', '<', 1}));
p.addParameter('BottomMargin', 0.07, @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'nonnegative', '<', 1}));
p.addParameter('WidthSpacing', 0.015, @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'nonnegative', '<', 1}));
p.addParameter('HeightSpacing', 0.035, @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'nonnegative', '<', 1}));

% Handle legacy positional arguments for backward compatibility
if ~isempty(varargin) && ~ischar(varargin{1}) && ~isstring(varargin{1})
    pos_args = cell(1, 6);
    for i = 1:min(length(varargin), 6)
        pos_args{i} = varargin{i};
    end
    
    if ~isempty(pos_args{1}), p.addParameter('LegacyWidth', ...
            pos_args{1}, @(x) true); end
    if ~isempty(pos_args{2}), p.addParameter('LegacyHeight', ...
            pos_args{2}, @(x) true); end
    if ~isempty(pos_args{3}), p.addParameter('LegacyLeftMargin', ...
            pos_args{3}, @(x) true); end
    if ~isempty(pos_args{4}), p.addParameter('LegacyBottomMargin', ...
            pos_args{4}, @(x) true); end
    if ~isempty(pos_args{5}), p.addParameter('LegacyWidthSpacing', ...
            pos_args{5}, @(x) true); end
    if ~isempty(pos_args{6}), p.addParameter('LegacyHeightSpacing', ...
            pos_args{6}, @(x) true); end
    
    p.parse(A);
    
    w = ifElse(~isempty(pos_args{1}), pos_args{1}, p.Results.Width);
    h = ifElse(~isempty(pos_args{2}), pos_args{2}, p.Results.Height);
    iwo = ifElse(~isempty(pos_args{3}), pos_args{3}, p.Results.LeftMargin);
    iho = ifElse(~isempty(pos_args{4}), pos_args{4}, p.Results.BottomMargin);
    iaws = ifElse(~isempty(pos_args{5}), pos_args{5}, p.Results.WidthSpacing);
    iahs = ifElse(~isempty(pos_args{6}), pos_args{6}, p.Results.HeightSpacing);
else
    p.parse(A, varargin{:});
    
    w = p.Results.Width;
    h = p.Results.Height;
    iwo = p.Results.LeftMargin;
    iho = p.Results.BottomMargin;
    iaws = p.Results.WidthSpacing;
    iahs = p.Results.HeightSpacing;
end

% Core layout calculations
mmw = @(x) x - (ceil(x/A(2))-1)*A(2);   % Column position from index
mmh = @(x) A(1) - ceil(x/A(2));          % Row position from index (inverted)

% Calculate individual axis dimensions
iaw = (w - (A(2)-1)*iaws) / A(2);        % Axis width
iah = (h - (A(1)-1)*iahs) / A(1);        % Axis height

% Create axes at calculated position
h = axes('Position', [iwo + (iaw + iaws)*(mmw(A(3)) - 1), ...
    iho + (iah + iahs)*mmh(A(3)), iaw, iah]);
end

function val = ifElse(condition, trueVal, falseVal)
    if condition
        val = trueVal;
    else
        val = falseVal;
    end
end