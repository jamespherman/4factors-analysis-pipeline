function metrics = calculate_waveform_metrics(waveform, fs)
% calculate_waveform_metrics - Computes key metrics from a spike waveform.
%
% This function analyzes a single mean waveform to extract features like
% the peak-to-trough duration.
%
% INPUTS:
%   waveform - A 1D vector representing the mean spike waveform.
%   fs       - The sampling rate of the waveform in Hz.
%
% OUTPUT:
%   metrics  - A struct containing the calculated waveform metrics:
%              - peak_trough_ms: The duration from the waveform's peak to
%                                its subsequent trough, in milliseconds.
%              - peak_val: The voltage value at the peak.
%              - trough_val: The voltage value at the trough.
%              - peak_idx: The sample index of the peak.
%              - trough_idx: The sample index of the trough.
%

% Ensure waveform is a row or column vector
if ~isvector(waveform)
    error('Input waveform must be a vector.');
end
waveform = waveform(:)'; % Ensure it's a row vector for consistency

% --- Find Peak and Trough ---
% Following the convention that the peak is the maximum positive
% deflection and the trough is the subsequent minimum deflection.

% Find the peak of the spike.
[peak_val, peak_idx] = max(waveform);

% Find the trough that occurs *after* the peak.
% We search from the peak index to the end of the waveform.
[trough_val, trough_idx_relative] = min(waveform(peak_idx:end));

% The index of the trough must be adjusted to be relative to the entire waveform.
trough_idx = peak_idx + trough_idx_relative - 1;

% --- Calculate Duration ---

% Duration in number of samples
duration_samples = trough_idx - peak_idx;

% Convert duration from samples to milliseconds
% duration_ms = (duration_in_samples / sampling_rate_in_Hz) * 1000
peak_trough_ms = (duration_samples / fs) * 1000;


% --- Store metrics in output struct ---
metrics.peak_trough_ms = peak_trough_ms;
metrics.peak_val = peak_val;
metrics.trough_val = trough_val;
metrics.peak_idx = peak_idx;
metrics.trough_idx = trough_idx;

end
