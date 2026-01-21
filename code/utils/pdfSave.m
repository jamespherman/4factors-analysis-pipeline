function pdfSave(fileName, paperSize, varargin)
%
% pdfSave(fileName, paperSize, [fig handle])
%
% Assuming you have a figure handle "fh" the "WYSIWYG" way to generate a
% PDF is:
%
% pdfSave(fileName, fh.Position(3:4)/72, fh)
%
% Uses exportgraphics with rasterization for smaller file sizes.
% Resolution is set to 150 DPI for good quality at reasonable size.

if nargin > 2
    fh = varargin{1};
else
    fh = gcf;
end

set(fh, 'PaperUnits', 'Inches', 'PaperSize', paperSize);
set(fh, 'PaperUnits', 'Normalized', 'PaperPosition', [0 0 1 1]);

% Use exportgraphics with rasterization for smaller file sizes
% 'ContentType', 'image' rasterizes the figure at specified resolution
exportgraphics(fh, fileName, 'Resolution', 150, 'ContentType', 'image');

end