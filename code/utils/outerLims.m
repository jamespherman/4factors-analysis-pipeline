function [xlims, ylims] = outerLims(ax)

flatt = @(x)x(:);

tempYlims = flatt(cell2mat(get(ax, 'YLim')));
tempXlims = flatt(cell2mat(get(ax, 'XLim')));

xlims = [min(tempXlims), max(tempXlims)];
ylims = [min(tempYlims), max(tempYlims)];

end