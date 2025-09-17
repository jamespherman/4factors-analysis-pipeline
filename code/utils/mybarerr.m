function [b, e] = mybarerr(x, y, yCI, varargin)
% [b, e] = mybarerr(x, y, yCI, Name, Value)
%
% Creates a bar plot with specified error bars. This function is a
% modern implementation that leverages MATLAB's built-in plotting
% functions and provides easy manipulation of graphics object appearances.
%
% INPUTS:
%   x                   - A vector specifying the x-coordinates of the bars.
%   y                   - A vector specifying the height of the bars.
%   yCI                 - A 2-column matrix specifying the confidence
%                         intervals for each bar. The first column holds
%                         the lower bounds, and the second column holds
%                         the upper bounds.
%
% OPTIONAL NAME-VALUE PAIR ARGUMENTS:
%   'BarWidth'          - Width of the bars. Default is 0.8.
%   'FaceColor'         - Color of the bars' faces. Can be a single color
%                         or an N-by-3 matrix for individual bar colors.
%                         Default is MATLAB's standard blue.
%   'EdgeColor'         - Color of the bars' edges. Default is 'none'.
%   'Orientation'       - Orientation, 'vertical' (default) or 'horizontal'.
%   'ErrorBarColor'     - Color of the error bars. Default is black.
%   'ErrorBarWidth'     - Line width of the error bars. Default is 1.5.
%   'ErrorBarCapSize'   - Size of the cap on the error bars. Default is 0.
%
% OUTPUTS:
%   b                   - The bar series graphics object.
%   e                   - The error bar series graphics object.
%
% EXAMPLE:
%   x = 1:5;
%   y = [20 35 40 55 70];
%   err = [5 8 10 6 9];
%   yCI = [y - err; y + err]';
%   figure;
%   [b, e] = mybarerr(x, y, yCI, 'FaceColor', 'b');
%   b.FaceAlpha = 0.5;
%   e.LineStyle = '--';
%   title('Modern Bar Plot with Error Bars');

    % Check for the correct number of input arguments.
    if nargin < 3
        error('mybarerr:NotEnoughInputs', ...
            'This function requires at least three input arguments.');
    end

    % Use inputParser for robust name-value pair argument handling.
    p = inputParser;

    % Define validation functions for arguments.
    isColorSpec = @(x) (ischar(x) && isvector(x)) || ...
        (isnumeric(x) && isvector(x) && length(x) == 3) || ...
        (isnumeric(x) && ismatrix(x) && size(x, 2) == 3);
    isValidOrientation = @(x) any(validatestring(x, ...
        {'vertical', 'horizontal'}));

    % Define default values for optional parameters.
    defaultBarWidth = 0.8;
    % Default FaceColor is MATLAB's standard blue.
    defaultFaceColor = [0 0.4470 0.7410];
    defaultEdgeColor = 'none';
    defaultOrientation = 'vertical';
    % Default ErrorBarColor is black.
    defaultErrorBarColor = [0 0 0];
    defaultErrorBarWidth = 1.5;
    defaultErrorBarCapSize = 0;

    % Add the parameters to the input parser.
    addParameter(p, 'BarWidth', defaultBarWidth, @isnumeric);
    addParameter(p, 'FaceColor', defaultFaceColor, isColorSpec);
    addParameter(p, 'EdgeColor', defaultEdgeColor, isColorSpec);
    addParameter(p, 'Orientation', defaultOrientation, isValidOrientation);
    addParameter(p, 'ErrorBarColor', defaultErrorBarColor, isColorSpec);
    addParameter(p, 'ErrorBarWidth', defaultErrorBarWidth, @isnumeric);
    addParameter(p, 'ErrorBarCapSize', defaultErrorBarCapSize, @isnumeric);

    % Parse the provided inputs.
    parse(p, varargin{:});

    % Store the parsed results in local variables.
    barWidth = p.Results.BarWidth;
    faceColor = p.Results.FaceColor;
    edgeColor = p.Results.EdgeColor;
    orientation = p.Results.Orientation;
    errBarColor = p.Results.ErrorBarColor;
    errBarWidth = p.Results.ErrorBarWidth;
    errBarCapSize = p.Results.ErrorBarCapSize;

    % Preserve the current hold state of the axes.
    wasHeld = ishold;
    if ~wasHeld
        hold on;
    end

    % yCI must be a matrix with two columns: [lower, upper].
    if size(yCI, 2) ~= 2
        error('mybarerr:InvalidCI', ...
            'yCI must be a matrix with two columns.');
    end
    
    % Ensure y is a column vector for calculations.
    if isrow(y)
        y = y';
    end

    % Calculate positive and negative error lengths from y and yCI.
    y_neg = y - yCI(:,1);
    y_pos = yCI(:,2) - y;

    % Plot bars and error bars based on the specified orientation.
    if strcmpi(orientation, 'vertical')
        % Plot vertical bars.
        b = bar(x, y, barWidth, 'EdgeColor', edgeColor);

        % Add vertical error bars.
        e = errorbar(x, y, y_neg, y_pos, 'LineStyle', 'none');
    else
        % Plot horizontal bars.
        b = barh(x, y, barWidth, 'EdgeColor', edgeColor);

        % Add horizontal error bars.
        e = errorbar(x, y, y_neg, y_pos, 'horizontal', ...
            'LineStyle', 'none');
    end

    % Set error bar properties.
    e.Color = errBarColor;
    e.LineWidth = errBarWidth;
    e.CapSize = errBarCapSize;

    % Apply face color to bars. This handles single or multi-color.
    isPerBarColor = size(faceColor, 1) == length(y) && ...
                    size(faceColor, 2) == 3;
    if isPerBarColor
        b.FaceColor = 'flat';
        b.CData = faceColor;
    else
        b.FaceColor = faceColor;
    end

    % Restore the original hold state of the axes.
    if ~wasHeld
        hold off;
    end

end