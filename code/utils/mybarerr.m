function [p, l] = mybarerr(x, y, yCI, varargin)
% [p, l] = mybarerr(x, y, yCI, Name, Value, ...)
%
% A modern implementation to plot bars with confidence intervals.
%
% This function creates bar-like plots where the extent of the bar
% represents a confidence interval (CI) and a line indicates the
% mean or central value. It is designed for easy manipulation of
% the resulting graphics objects.
%
% EXAMPLE:
%   x = 1:5;
%   y = rand(1, 5) * 5;
%   err = rand(1, 5);
%   yCI = [y - err; y + err]';
%   figure;
%   mybarerr(x, y, yCI);
%
% INPUTS:
%   x       - A vector of x-positions for the bars.
%   y       - A vector of y-values (e.g., means).
%   yCI     - An n-by-2 matrix of confidence intervals, [lower upper].
%
% NAME-VALUE PAIRS:
%   'yCI2' - A second, inner CI to plot. Default is [].
%   'Colors' - An n-by-3 matrix of RGB colors. Default is generated
%              from the 'parula' colormap.
%   'BarWidth' - Scalar value to control the width of the bars.
%                Default is auto-calculated based on x-spacing.
%   'Orientation' - 'vertical' (default) or 'horizontal'.
%
% OUTPUTS:
%   p       - Handles to the patch objects (the bars). If yCI2 is
%             used, p is a struct with fields p.outer and p.inner.
%   l       - Handles to the line objects (the means).

% --- Input Parsing ---
parser = inputParser;

% Validation functions
isNumericReal = @(a) isnumeric(a) && isreal(a);
isValidColors = @(c) ismatrix(c) && size(c, 2) == 3;
isValidWidth = @(w) isscalar(w) && isNumericReal(w);

% Add required arguments
addRequired(parser, 'x', isNumericReal);
addRequired(parser, 'y', isNumericReal);
addRequired(parser, 'yCI', isNumericReal);

% Add optional name-value pair arguments
addParameter(parser, 'yCI2', [], isNumericReal);
addParameter(parser, 'Colors', [], isValidColors);
addParameter(parser, 'BarWidth', [], isValidWidth);
addParameter(parser, 'Orientation', 'vertical', @isValidOrientation);

% Execute the parser
parse(parser, x, y, yCI, varargin{:});

% --- Extract values from parser ---
args = parser.Results;
x = args.x(:);
y = args.y(:);
yCI = args.yCI;
yCI2 = args.yCI2;
colors = args.Colors;
barWidth = args.BarWidth;
orientation = lower(args.Orientation);

% --- Data and Color Setup ---
nBars = length(x);

% Make sure inputs are the correct size
if length(y) ~= nBars || size(yCI, 1) ~= nBars
    error('Inputs x, y, and yCI must have the same number of elements.');
end
if ~isempty(yCI2) && size(yCI2, 1) ~= nBars
    error('Input yCI2 must have the same number of elements as x.');
end

% If no colors are provided, generate them using a default colormap.
if isempty(colors)
    colors = parula(nBars);
end

% If only one color is given for many bars, repeat it for all bars.
if size(colors, 1) == 1 && nBars > 1
    colors = repmat(colors, nBars, 1);
end

% --- Bar Width Setup ---
if isempty(barWidth)
    if length(x) > 1
        % Default width is half the minimum distance between x-points.
        sorted_x = sort(x);
        halfWidth = 0.4 * min(diff(sorted_x));
    else
        % Default width for a single bar.
        halfWidth = 0.4;
    end
else
    % Use user-specified width, converting radius to half-width.
    halfWidth = barWidth / 2;
end

% --- Plotting ---

% Get current hold state and then turn hold on for plotting.
wasHeld = ishold;
hold on;

% Pre-allocate graphics handle arrays for efficiency.
p_outer = gobjects(nBars, 1);
if ~isempty(yCI2)
    p_inner = gobjects(nBars, 1);
end
l = gobjects(nBars, 1);

% Loop through each bar to create the graphics objects individually.
for i = 1:nBars
    % Define colors for the elements of the current bar.
    outerColor = colors(i, :);
    innerColor = outerColor * 0.7; % Darker for the inner bar
    lineColor = outerColor * 0.5;  % Darkest for the mean line

    % Define bar vertices and line coordinates based on orientation.
    if strcmp(orientation, 'vertical')
        X1 = [x(i)-halfWidth, x(i)+halfWidth, x(i)+halfWidth, x(i)-halfWidth];
        Y1 = [yCI(i,1), yCI(i,1), yCI(i,2), yCI(i,2)];
        X_line = [x(i) - 0.9*halfWidth, x(i) + 0.9*halfWidth];
        Y_line = [y(i), y(i)];
    else % horizontal
        Y1 = [x(i)-halfWidth, x(i)+halfWidth, x(i)+halfWidth, x(i)-halfWidth];
        X1 = [yCI(i,1), yCI(i,1), yCI(i,2), yCI(i,2)];
        Y_line = [x(i) - 0.9*halfWidth, x(i) + 0.9*halfWidth];
        X_line = [y(i), y(i)];
    end

    % Draw the outer bar for the primary confidence interval.
    p_outer(i) = patch(X1, Y1, outerColor, 'EdgeColor', 'none');

    % If a second CI is provided, draw it on top of the outer bar.
    if ~isempty(yCI2)
        if strcmp(orientation, 'vertical')
            Y2 = [yCI2(i,1), yCI2(i,1), yCI2(i,2), yCI2(i,2)];
            p_inner(i) = patch(X1, Y2, innerColor, 'EdgeColor', 'none');
        else % horizontal
            X2 = [yCI2(i,1), yCI2(i,1), yCI2(i,2), yCI2(i,2)];
            p_inner(i) = patch(X2, Y1, innerColor, 'EdgeColor', 'none');
        end
    end

    % Draw the mean line, which is slightly indented from the bar edges.
    l(i) = line(X_line, Y_line, 'Color', lineColor, 'LineWidth', 3);
end

% --- Finalize ---

% Restore the original hold state of the plot.
if ~wasHeld
    hold off;
end

% Set the output patch handle. If a second CI was used, return a
% struct containing handles for both inner and outer patches.
if ~isempty(yCI2)
    p.outer = p_outer;
    p.inner = p_inner;
else
    p = p_outer;
end
end

function TF = isValidOrientation(s)
    % Local function to validate the Orientation argument.
    isText = ischar(s) || isstring(s);
    isValid = any(strcmpi(s, {'vertical', 'horizontal'}));
    TF = isText && isValid;
end
