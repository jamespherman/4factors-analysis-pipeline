function subregion = parse_snc_subregion(grid_hole_str, brain_area)
% PARSE_SNC_SUBREGION Determine SNc subregion from grid coordinates.
%
% This function parses the grid hole coordinate string from the session
% manifest and determines the SNc subregion based on the Y-coordinate.
%
% INPUT:
%   grid_hole_str: String in format '(+X, +Y)' from manifest
%                  (e.g., '(+5, +3)', '(+4,+0.5)')
%   brain_area:    String indicating brain area ('SC' or 'SNc')
%
% OUTPUT:
%   subregion: 'rvmSNc', 'cdlSNc', or '' (empty for non-SNc sessions)
%
% CLASSIFICATION RULE (only applied when brain_area == 'SNc'):
%   Y >= 3 -> 'rvmSNc' (rostral ventromedial)
%   Y < 3  -> 'cdlSNc' (caudal dorsolateral)
%
% Author: Claude Code
% Date: 2026-01-08

% Return empty string for non-SNc sessions
if ~strcmp(brain_area, 'SNc')
    subregion = '';
    return;
end

% Parse coordinates using sscanf which handles format variations
% Supports: '(+5, +3)', '(+4,+0.5)', '(-2,+1.5)', '(+3.5, +2.5)'
coords = sscanf(grid_hole_str, '(%f, %f)');

% Handle format without space after comma
if length(coords) < 2
    coords = sscanf(grid_hole_str, '(%f,%f)');
end

if length(coords) >= 2
    grid_y = coords(2);
    if grid_y >= 3
        subregion = 'rvmSNc';
    else
        subregion = 'cdlSNc';
    end
else
    subregion = 'unknown';
    warning('parse_snc_subregion:parseFailed', ...
        'Could not parse grid_hole: %s', grid_hole_str);
end

end
