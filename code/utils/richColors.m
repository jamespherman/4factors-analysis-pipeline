function colors = richColors

% define "main" colors
rgbValues = [...
     255 194 10;  % Bright Yellow
     12 123 220;  % Vivid Blue
     153 79 0;    % Dark Orange
     0 108 209;   % Deep Sky Blue
     26 255 26;   % Lime Green
     75 0 146;    % Deep Purple
     254 254 98;  % Light Yellow
     211 95 183;  % Medium Pink
     225 190 106; % Light Orange
     64 176 166;  % Teal
     0 90 181;    % Royal Blue
     220 50 32;   % Bright Red
     230 97 90;   % Coral
     93 58 155;   % Medium Purple
     26 133 255;  % Bright Blue
     212 17 89;   % Deep Pink
    ];

% define "my" order
myOrder = [6; 14; 11; 4; 2; 15; 10; 5; 7; 9; 3; 1; 13; 12; 16; 8];

% reorder and return
colors = rgbValues(myOrder, :) / 255;