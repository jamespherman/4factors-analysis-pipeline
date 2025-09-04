function h = mySubPlot(A, varargin)
% MYSUBPLOT Create custom subplot with precise positioning control
%
%   h = mySubPlot(A [, 'PropertyName', PropertyValue, ...])
%
% A is a 1x3 or 3x1 vector specifying the subplot parameters:
%   A(1): Number of rows in the subplot grid
%   A(2): Number of columns in the subplot grid
%   A(3): Index of the current subplot (fills row-wise)
%
% Optional parameters (name-value pairs):
%   'Width'      - Total width of subplot grid (default: 0.85)
%   'Height'     - Total height of subplot grid (default: 0.85)
%   'LeftMargin' - Left margin before first column (default: 0.07)
%   'BottomMargin' - Bottom margin before first row (default: 0.07)
%   'WidthSpacing' - Width gap between subplots (default: 0.015)
%   'HeightSpacing' - Height gap between subplots (default: 0.035)
%
% Example:
%   mySubPlot([2,3,1]) % Upper-left subplot in a 2x3 grid
%   mySubPlot([2,3,4], 'Width', 0.9, 'HeightSpacing', 0.05)
%
% Original by jph - 1/23/2015, modernized with inputParser - 2/25/2025
% Create input parser
p = inputParser;
p.FunctionName = 'mySubPlot';
p.CaseSensitive = false;
p.KeepUnmatched = true;
% Validate subplot grid specification
validateA = @(x) validateattributes(x, {'numeric'}, ...
    {'vector', 'numel', 3, 'positive', 'integer'});
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
% For backward compatibility, also accept positional arguments
if ~isempty(varargin) && ~ischar(varargin{1}) && ...
        ~isstring(varargin{1})
    % Handle legacy positional arguments
    pos_args = cell(1, 6);
    for i = 1:min(length(varargin), 6)
        pos_args{i} = varargin{i};
    end
    
    % Parse positional arguments if provided
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
    
    % Extract parameters, prioritizing legacy positional arguments
    w = ifElse(~isempty(pos_args{1}), pos_args{1}, p.Results.Width);
    h = ifElse(~isempty(pos_args{2}), pos_args{2}, p.Results.Height);
    iwo = ifElse(~isempty(pos_args{3}), pos_args{3}, ...
        p.Results.LeftMargin);
    iho = ifElse(~isempty(pos_args{4}), pos_args{4}, ...
        p.Results.BottomMargin);
    iaws = ifElse(~isempty(pos_args{5}), pos_args{5}, ...
        p.Results.WidthSpacing);
    iahs = ifElse(~isempty(pos_args{6}), pos_args{6}, ...
        p.Results.HeightSpacing);
else
    % Modern name-value pair parsing
    p.parse(A, varargin{:});
    
    % Extract parameters
    w = p.Results.Width;
    h = p.Results.Height;
    iwo = p.Results.LeftMargin;
    iho = p.Results.BottomMargin;
    iaws = p.Results.WidthSpacing;
    iahs = p.Results.HeightSpacing;
end

% In-line function definitions (core layout algorithm)
mmw = @(x) x - (ceil(x/A(2))-1)*A(2);
mmh = @(x) A(1) - ceil(x/A(2));

% Calculate the individual axis width and height based on number of axes
% and specified geometric values.
iaw = (w - (A(2)-1)*iaws) / A(2);
iah = (h - (A(1)-1)*iahs) / A(1);

% Make axes
h = axes('Position', [iwo + (iaw + iaws)*(mmw(A(3)) - 1), ...
    iho + (iah + iahs)*mmh(A(3)), iaw, iah]);
end

% Helper function to handle conditional assignment
function val = ifElse(condition, trueVal, falseVal)
    if condition
        val = trueVal;
    else
        val = falseVal;
    end
end