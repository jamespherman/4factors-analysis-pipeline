function pngSave(fileName, paperSize, varargin)
%
% pngSave(fileName, paperSize, [fig handle])
%
% 

if nargin > 2
    fh = varargin{1};
else
    fh = gcf;
end

set(fh, 'PaperUnits', 'Inches', 'PaperSize', paperSize);
set(fh, 'PaperUnits', 'Normalized', 'PaperPosition', [0 0 1 1]);
saveas(fh, fileName, 'png')

end