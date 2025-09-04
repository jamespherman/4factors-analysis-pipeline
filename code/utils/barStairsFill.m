function hOut = barStairsFill(x, y1, y2)

% hOut = barStairsFill(x, y1, y2)
%
% plot two traces (y1, and y2) in my "stairs" style, and fill in the area
% between the two traces with gray; return graphics handles to the fill
% object, hOut(1), and the two stairs traces (hOut(2:3);

xDiff = mean(diff(x))/2;
xStairs     = [x(:)'-xDiff, x(end)+xDiff];
yStairs1    = [y1(:)', y1(end)];
yStairs2    = [y2(:)', y2(end)];
[xx, yy1]   = stairs(xStairs, yStairs1);
[~, yy2]    = stairs(xStairs, yStairs2);

yAll    = [yy1(:)'; yy2(:)'];
xFill   = [xx(:)', fliplr(xx(:)'), xx(1)];
yMin    = min(yAll);
yMax    = max(yAll);
yFill   = [yMin, fliplr(yMax), yMin(1)];

hold on
hOut(1) = fill(xFill, yFill, 0.7 * [1 1 1], 'EdgeColor', 'none');
hOut(2) = barStairs(x, y1, false);
hOut(3) = barStairs(x, y2, false);

end