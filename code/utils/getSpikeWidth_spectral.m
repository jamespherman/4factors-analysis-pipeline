function duration_ms = getSpikeWidth_spectral(waveform, fs)
% Calculates spike width based on the inverse of the dominant frequency
% in the waveform's power spectrum.
%
% INPUTS:
%   waveform - A 1D vector containing the extracellular spike waveform.
%   fs       - The sampling rate of the recording in Hz.
%
% OUTPUT:
%   duration_ms - The calculated duration of the spike in milliseconds.

% Ensure the waveform is a column vector for consistency
if isrow(waveform)
    waveform = waveform';
end

% --- Step 1: Detrend and prepare the signal ---
% It's good practice to remove any DC offset before the FFT.
waveform = detrend(waveform, 'constant'); % Removes the mean value

% --- Step 2: Calculate the power spectrum using FFT ---
L = length(waveform);       % Length of signal
Y = fft(waveform);          % Compute the fast fourier transform
P2 = abs(Y/L);              % Two-sided spectrum
P1 = P2(1:floor(L/2)+1);    % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);

% --- Step 3: Find the dominant frequency ---
% Create the frequency vector for the x-axis of the spectrum
f = fs*(0:(L/2))/L;

% Find the maximum power, IGNORING the DC component (the first element)
% The DC component at f=0 is often the largest but isn't relevant for
% the waveform's oscillatory shape.
[~, max_idx] = max(P1(2:end)); % Search from the second element onwards

% Get the frequency corresponding to the max power
% We add 1 to the index because we started our search at P1(2)
dominant_frequency = f(max_idx + 1);

% --- Step 4: Calculate duration in seconds ---
duration_sec = 1 / dominant_frequency;

% --- Step 5: Convert to milliseconds ---
duration_ms = duration_sec * 1000;

end