function h = barStairs(x, y, varargin)

% h = barStairs(x, y, [wantBar])
%
% h(1) is the bar series graphics object handle
% h(2) is the stairs graphics object handle
%

if nargin > 2
    wantBar = varargin{1};
else
    wantBar = true;
end

xdiff = mean(diff(x))/2;

xStairs = [x(:)'-xdiff, x(end)+xdiff];
yStairs = [y(:)', y(end)];

hold on;
hInd = 1;
if wantBar
    h(hInd) = bar(x,y,1,'EdgeColor', 'None', 'FaceColor', 0.75*[1 1 1]);
    hInd = hInd + 1;
end
h(hInd) = stairs(xStairs, yStairs, 'Color', [0 0 0], 'LineWidth', 2);

end