%% richColors.m
%
% Provides a standardized palette of 16 perceptually distinct colors for
% use in plots and figures.
%
% This function returns a predefined set of colors that are visually
% distinguishable, making them suitable for categorical data representation.
% The color palette includes a range of hues, such as purples, blues,
% greens, yellows, oranges, and reds.
%
% The order of colors is arranged to maximize contrast between adjacent
% colors in the list. For instance, 'Deep Purple' is followed by 'Medium
% Purple', but then transitions to 'Royal Blue' to ensure a clear visual
% separation. When selecting colors for a plot, consider the desired level
% of contrast. For high contrast, select colors from distant positions in
% the list. For a more harmonious or related color scheme, select colors
% that are adjacent or have similar names (e.g., 'Vivid Blue' and 'Bright
% Blue'). The RGB values can also be inspected to gauge similarity.
%
%
% Usage:
%   colors = richColors()
%   colors = richColors('matrix')
%   colors = richColors('table')
%
% Inputs:
%   output_format (string, optional): Specifies the format of the output.
%       - 'matrix' (default): Returns an Nx3 matrix of normalized RGB
%                             triplets, where N is the number of colors.
%       - 'table': Returns a table with columns for 'Name' and 'RGB'
%                  values (normalized).
%
% Outputs:
%   colors: An Nx3 matrix or a table containing the color palette,
%           depending on the specified output_format.
%
% Author: Your Name
% Date: 2025-09-26
%

function colors = richColors(output_format)
    % Define the color data as a cell array. Each row contains a color
    % name and its corresponding [R, G, B] values in the range [0, 255].
    color_data = { ...
        'Deep Purple',      [75 0 146]; ...
        'Medium Purple',    [93 58 155]; ...
        'Royal Blue',       [0 90 181]; ...
        'Deep Sky Blue',    [0 108 209]; ...
        'Vivid Blue',       [12 123 220]; ...
        'Bright Blue',      [26 133 255]; ...
        'Teal',             [64 176 166]; ...
        'Lime Green',       [26 255 26]; ...
        'Light Yellow',     [254 254 98]; ...
        'Bright Yellow',    [255 194 10]; ...
        'Light Orange',     [225 190 106]; ...
        'Dark Orange',      [153 79 0]; ...
        'Bright Red',       [220 50 32]; ...
        'Coral',            [230 97 90]; ...
        'Deep Pink',        [212 17 89]; ...
        'Medium Pink',      [211 95 183] ...
    };

    % Default to 'matrix' if no format is specified
    if nargin < 1
        output_format = 'matrix';
    end

    % Return the requested format
    switch output_format
        case 'table'
            color_table = cell2table(color_data, ...
                'VariableNames', {'Name', 'RGB'});
            % Normalize RGB values
            color_table.RGB = color_table.RGB / 255;
            colors = color_table;
        case 'matrix'
            % Return just the RGB values, normalized
            colors = cell2mat(color_data(:,2)) / 255;
        otherwise
            error("Invalid output_format. Use 'matrix' or 'table'.");
    end
end