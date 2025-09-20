% --- Function Definition ---
function offset_array = offsetArrayWaveform(waveforms, percentile_margin)
% CREATE_OFFSET_ARRAY Offsets multi-channel data for clean plotting.
%
%   offset_array = CREATE_OFFSET_ARRAY(waveforms, percentile_margin)
%   takes an nChannels x nSamples array `waveforms` and calculates a
%   vertical offset to separate the channels for plotting. The offset
%   is determined by the maximum range between specified percentiles.
%
%   INPUTS:
%   waveforms         - An [nChannels x nSamples] numerical array of
%                       voltage data.
%   percentile_margin - A scalar proportion between 0 and 0.5. For
%                       example, 0.01 will use the 1st and 99th
%                       percentiles to define the range.
%
%   OUTPUT:
%   offset_array      - A [nChannels x nSamples] array where each
%                       channel (row) has been vertically offset from
%                       the one below it.

% --- 1. Input Validation ---
if ~isnumeric(waveforms) || ~ismatrix(waveforms)
    error('Input `waveforms` must be a 2D numerical array.');
end
if ~isscalar(percentile_margin) || percentile_margin < 0 || ...
        percentile_margin >= 0.5
    error('Input `percentile_margin` must be a scalar between 0 and 0.5.');
end

% --- 2. Calculate Offset Value ---
nChannels= size(waveforms, 1);

% Define the lower and upper percentile bounds based on the margin
prctile_bounds = [percentile_margin, 1 - percentile_margin] * 100;

% Calculate the percentile values for each channel along the samples 
% dimension (dim 2). This returns an nChannels x 2 matrix.
channel_percentiles = prctile(waveforms, prctile_bounds, 2);

% Calculate the voltage range for each channel based on these percentiles
channel_ranges = channel_percentiles(:, 2) - channel_percentiles(:, 1);

% The final offset value is the maximum of these ranges
offset_value = max(channel_ranges);

% --- 3. Apply Offset ---
% Create a column vector of cumulative offsets, one for each channel.
% We create it descending so channel 1 is at the top of the plot.
offsets = (nChannels-1:-1:0)' * offset_value;

% Use implicit expansion (broadcasting) to add the offset vector to
% each sample in the corresponding channel. This is more efficient
% than repmat for modern MATLAB versions.
offset_array = waveforms + offsets;

end