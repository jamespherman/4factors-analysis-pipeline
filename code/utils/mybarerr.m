function [fo, po, io] = mybarerr(x, y, yCI, varargin)
%MYBARERR  Draw box-style bars whose extent is a confidence interval.
%
%   [fo, po]     = MYBARERR(x, y, yCI)
%   [fo, po]     = MYBARERR(x, y, yCI, yi)
%   [fo, po]     = MYBARERR(x, y, yCI, yi, color)
%   [fo, po]     = MYBARERR(x, y, yCI, yi, color, radius)
%   [fo, po]     = MYBARERR(x, y, yCI, yi, color, radius, horiz)
%   [fo, po, io] = MYBARERR(...)
%
%   MYBARERR draws one filled rectangle ("bar") per element of x.
%   Unlike a normal bar, each bar spans a confidence interval (yCI)
%   rather than rising from zero, and a thick line marks the central
%   value (y). It is a compact way to show mean/median +/- CI. Call
%   it once for a whole group, or once per bar inside a loop when you
%   need per-bar colors or widths.
%
% REQUIRED INPUTS
%   x    1xN or Nx1. Category-axis positions of the bar centers.
%   y    1xN or Nx1. Central value per bar (drawn as a thick line).
%   yCI  Nx2 (or 2xN, auto-transposed). [lower upper] extent of each
%        bar, one row per bar. These are ABSOLUTE axis values, not
%        +/- offsets from y.
%
% OPTIONAL INPUTS (positional; pass [] to skip and keep the default)
%   yi     Second, inner interval. Its role depends on its shape:
%            * Nx2 -> inner CI; a nested box is drawn inside yCI
%              (e.g. a mean CI inside a 95% spread). Enables 2-box
%              mode and returns inner handles in io.
%            * scalar in [0,1] -> whitening level for the auto
%              rainbow palette (see color): 0 = vivid, 1 = white.
%            * [] -> no inner box and no rainbow whitening. (default)
%   color  Bar fill color. Accepts:
%            * 1x3 RGB     -> the same fill for every bar.
%            * Nx3 RGB     -> one fill per bar; each center line is a
%              darkened version of that bar's fill.
%            * scalar true -> auto rainbow palette (hsv), whitened by
%              yi when yi is scalar.
%            * [] or 0     -> default medium gray. (default)
%   radius Half-width of each bar on the category axis. [] -> auto:
%          40% of the smallest gap between x values, capped at 0.4;
%          a single bar uses 0.4. (default [])
%   horiz  false/0 -> vertical bars, value on the y-axis. (default)
%          true/1  -> horizontal bars, value on the x-axis.
%
% OUTPUTS
%   fo   1xN graphics handles to the (outer) fill objects.
%   po   1xN graphics handles to the central-value lines.
%   io   1xN handles to the inner fills in 2-box mode; empty
%        placeholders otherwise.
%
% BEHAVIOR NOTES
%   * Drawing is additive. MYBARERR turns HOLD ON if needed and
%     restores the prior hold state on exit, so set axis limits after
%     calling. It never creates a figure or clears the axes.
%   * Each fill's EdgeColor matches its face, so bars look
%     borderless.
%   * yCI must have one row per bar; a 2xN yCI is auto-transposed.
%
% EXAMPLES
%   % One gray bar: mean mu with 95% CI [lo hi] at x = 1.
%   mybarerr(1, mu, [lo hi]);
%
%   % Three bars at x = 1:3, each its own RGB fill, width 0.3.
%   for k = 1:3
%       mybarerr(k, mu(k), ci(k,:), [], rgb(k,:), 0.3);
%   end
%
%   % Auto rainbow, lightly whitened (0.2), five bars.
%   mybarerr(1:5, mu, ci, 0.2, true);
%
%   % Nested boxes: 95% spread (yCI) with a mean CI (yi) inside.
%   mybarerr(1:3, mu, spread95, meanCI);
%
%   % Horizontal bars (value on the x-axis).
%   mybarerr(1:3, mu, ci, [], [], [], true);
%
% LLM / AGENT USAGE NOTES
%   * Positional arguments only; there are no name-value pairs.
%   * yCI holds absolute [lower upper] bounds, NOT +/- half-widths.
%     For half-widths h, pass [y-h, y+h].
%   * The 4th arg (yi) is overloaded: an Nx2 matrix is an inner
%     interval, a scalar is a rainbow whitening level. Do not pass a
%     scalar expecting it to be an interval.
%   * For per-bar colors, pass an Nx3 matrix as `color`; you do NOT
%     also set the scalar true flag in that case.
%   * To give each bar a distinct color AND a manual width, call in a
%     loop with one bar per call, as in the per-bar RGB example.
%
%   See also FILL, PLOT, PATCH, ERRORBAR, HSV.

% Divisor used to darken a fill color into its center-line color.
lusv = 1.5;

% Default fill (medium gray) and its darkened center-line color.
defFill = [0.6 0.6 0.6];
defLine = defFill / lusv;

% Number of bars to draw.
nBoxes = numel(x);

% Unpack optional positional inputs with their defaults.
yi = []; colorArg = []; brm = []; horiz = 0;
if numel(varargin) >= 1, yi       = varargin{1}; end
if numel(varargin) >= 2, colorArg = varargin{2}; end
if numel(varargin) >= 3, brm      = varargin{3}; end
if numel(varargin) >= 4, horiz    = varargin{4}; end

% Ensure yCI has one row per bar (accept a transposed 2xN input).
if size(yCI, 1) ~= nBoxes
    yCI = yCI';
end

% 2-box mode is on only when yi is a (non-scalar) interval matrix.
inci = ~isempty(yi) && ~isscalar(yi);

% Match the inner interval's orientation to yCI when in 2-box mode.
if inci && size(yi, 1) ~= nBoxes
    yi = yi';
end

% Resolve one fill color and one line color per bar.
fillColors = repmat(defFill, nBoxes, 1);
lineColors = repmat(defLine, nBoxes, 1);
if ~isempty(colorArg)
    if isscalar(colorArg)
        % A truthy scalar requests the auto rainbow palette, but only
        % a scalar yi supplies a whitening level; otherwise stay gray.
        if colorArg && isscalar(yi)
            [fillColors, lineColors] = rainbowColors(nBoxes, yi, lusv);
        end
    elseif size(colorArg, 1) == 1
        % A single RGB triplet fills every bar; darken it for lines.
        fillColors = repmat(colorArg, nBoxes, 1);
        lineColors = repmat(colorArg / lusv, nBoxes, 1);
    else
        % One RGB triplet per bar; darken each for its center line.
        fillColors = colorArg;
        lineColors = colorArg / lusv;
    end
end

% Choose the bar half-width: manual, auto from spacing, or a fixed
% fallback for a lone bar.
if nBoxes > 1 || ~isempty(brm)
    if isempty(brm)
        br = min(0.4 * min(diff(sort(x(:)))), 0.4);
    else
        br = brm;
    end
else
    br = 0.4;
end

% Draw additively; remember the hold state so we can restore it.
wasHeld = ishold;
if ~wasHeld
    hold on
end

% Preallocate handle arrays (io stays empty unless in 2-box mode).
fo = gobjects(1, nBoxes);
po = gobjects(1, nBoxes);
io = gobjects(1, nBoxes);

for i = 1:nBoxes
    fc = fillColors(i, :);
    lc = lineColors(i, :);
    if inci
        [fo(i), po(i), io(i)] = drawNested(x(i), y(i), ...
            yCI(i, :), yi(i, :), br, horiz, fc, lc, lusv);
    else
        [fo(i), po(i)] = drawSingle(x(i), y(i), yCI(i, :), ...
            br, horiz, fc, lc);
    end
end

% Restore the caller's original hold state.
if ~wasHeld
    hold off
end
end


function [fo, po] = drawSingle(xc, yv, ci, br, horiz, fc, lc)
%DRAWSINGLE  One CI box with a central-value line.

% Left/right edges of the box on the category axis.
x1 = xc - br;
x2 = xc + br;

% Box corners: category-axis span and value-axis [lo hi hi lo].
boxCat = [x2 x2 x1 x1];
boxVal = [ci(1) ci(2) ci(2) ci(1)];

% Center line, inset slightly from the box edges.
lineCat = [x1 + br/10, x2 - br/10];
lineVal = [yv yv];

% Swap axes for horizontal bars.
if ~horiz
    fo = fill(boxCat, boxVal, fc, 'EdgeColor', fc);
    po = plot(lineCat, lineVal, 'LineWidth', 3, 'Color', lc);
else
    fo = fill(boxVal, boxCat, fc, 'EdgeColor', fc);
    po = plot(lineVal, lineCat, 'LineWidth', 3, 'Color', lc);
end
end


function [fo, po, io] = drawNested(xc, yv, ci, inr, br, horiz, ...
    fc, lc, lusv)
%DRAWNESTED  Outer CI box, inner CI box, and a central-value line.

% Outer box edges and a narrower inner box (55% of the half-width).
x1 = xc - br;
x2 = xc + br;
ir = 0.55 * br;
xi1 = xc - ir;
xi2 = xc + ir;

% Outer and inner box corners on category and value axes.
outerCat = [x2 x2 x1 x1];
outerVal = [ci(1) ci(2) ci(2) ci(1)];
innerCat = [xi2 xi2 xi1 xi1];
innerVal = [inr(1) inr(2) inr(2) inr(1)];

% Center line, inset slightly from the outer box edges.
lineCat = [x1 + br/10, x2 - br/10];
lineVal = [yv yv];

% Inner box is a darkened shade of the outer fill.
inColor = fc / lusv;

% Swap axes for horizontal bars.
if ~horiz
    fo = fill(outerCat, outerVal, fc, 'EdgeColor', fc);
    io = fill(innerCat, innerVal, inColor, 'EdgeColor', inColor);
    po = plot(lineCat, lineVal, 'LineWidth', 3, 'Color', lc);
else
    fo = fill(outerVal, outerCat, fc, 'EdgeColor', fc);
    io = fill(innerVal, innerCat, inColor, 'EdgeColor', inColor);
    po = plot(lineVal, lineCat, 'LineWidth', 3, 'Color', lc);
end
end


function [fc, lc] = rainbowColors(n, t, lusv)
%RAINBOWCOLORS  Whitened hsv palette; t in [0,1], 0 vivid, 1 white.

% t < 1 whitens by t; t >= 1 falls back to a light pastel.
if t < 1
    w = 1 - t;
    a = t;
else
    w = 0.2;
    a = 0.8;
end

% Blend the hsv wheel toward white, then darken for the lines.
fc = hsv(n) * w + a;
lc = fc / lusv;
end
