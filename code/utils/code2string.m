function out = code2string(values, times, codes)
% CODE2STRING Convert numeric trial codes to their string representations
%   out = code2string(values, times, codes) converts numeric codes in 'values'
%   to their corresponding field names from the 'codes' struct. Also includes
%   timing information if provided.
%
% Inputs:
%   values - Array of numeric codes to convert to strings
%   times  - Optional timing values associated with each code (can be empty)
%   codes  - Struct where field names are string labels and values are codes
%
% Output:
%   out    - nx2 cell array where:
%            Column 1: String labels matching each input code
%            Column 2: Timing values (if provided) or NaN

% Get all field names from the codes struct
fieldNames = fieldnames(codes);

% Extract all code values from struct using structfun
% This creates a vector of all the numeric codes
codeVals = structfun(@(x) x, codes);

% Initialize output cell array
nValues = length(values);
out = cell(nValues, 2);

% Convert each numeric code to its string representation
for i = 1:nValues
    if ~any(codeVals == values(i))
        % If code isn't found in the struct, mark as 'not found'
        outString = 'not found';
    else
        % Get the field name corresponding to this code value
        outString = fieldNames{codeVals == values(i)};
    end
    out{i, 1} = outString;
end

% Fill in timing information
if ~isempty(times)
    % If timing values provided, include them in second column
    out(:,2) = num2cell(times);
else
    % Otherwise fill with NaN
    out(:,2) = {NaN};
end
end