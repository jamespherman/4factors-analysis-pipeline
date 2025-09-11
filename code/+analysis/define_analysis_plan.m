function analysisPlan = define_analysis_plan()
% DEFINE_ANALYSIS_PLAN creates a struct with parameters for the 4-factors analysis.
%
%   This function centralizes the definition of the analysis plan,
%   making it easy to modify analysis parameters in one place.
%
%   Usage:
%   plan = analysis.define_analysis_plan()
%
%   Outputs:
%   analysisPlan - A struct with the following fields:
%       .alignEvents: A cell array of event names to align spike data to.
%       .timeWindow:  A 1x2 vector specifying the [min, max] time in
%                     seconds around the alignment event.
%       .binWidth:    A scalar specifying the bin size in seconds for
%                     creating histograms.

% Define the events of interest for alignment
analysisPlan.alignEvents = {'targetOn', 'saccadeOnset', 'reward'};

% Define the time window for analysis around each event
% (e.g., -1.0 to +2.0 seconds)
analysisPlan.timeWindow = [-1.0, 2.0];

% Define the width of the bins for spike rate calculation
% (e.g., 50 milliseconds)
analysisPlan.binWidth = 0.05;

end
