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
% Anchor on the largest absolute deflection from baseline. This handles
% both standard (trough-first) and inverted (peak-first) spikes.

% Find the index of the largest absolute deflection
[~, idx_abs_max] = max(abs(waveform));

% Check if this primary deflection is negative (trough) or positive (peak)
if waveform(idx_abs_max) < 0
    % The primary deflection is a trough.
    trough_idx = idx_abs_max;
    trough_val = waveform(trough_idx);

    % Find the subsequent positive peak.
    [peak_val, temp_idx] = max(waveform(trough_idx:end));
    peak_idx = temp_idx + trough_idx - 1; % Correct for sub-array indexing

else
    % The primary deflection is a peak.
    peak_idx = idx_abs_max;
    peak_val = waveform(peak_idx);

    % Find the subsequent negative trough.
    [trough_val, temp_idx] = min(waveform(peak_idx:end));
    trough_idx = temp_idx + peak_idx - 1; % Correct for sub-array indexing
end

% --- Calculate Duration ---

% Convert duration from samples to milliseconds
% duration_ms = (duration_in_samples / sampling_rate_in_Hz) * 1000
peak_trough_ms = getSpikeWidth_spectral(waveform, fs);

% --- Store metrics in output struct ---
metrics.peak_trough_ms = peak_trough_ms;
metrics.peak_val = peak_val;
metrics.trough_val = trough_val;
metrics.peak_idx = peak_idx;
metrics.trough_idx = trough_idx;

end
