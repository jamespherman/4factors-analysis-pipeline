function [a, ci, tpfp, sig] = arrayROC(x1, x2, varargin)
% arrayROC  Compute ROC area under curve for paired array columns
%
%   [a, ci, tpfp, sig] = arrayROC(x1, x2)
%   [a, ci, tpfp, sig] = arrayROC(x1, x2, reps)
%   [a, ci, tpfp, sig] = arrayROC(x1, x2, reps, alpha)
%   [a, ci, tpfp, sig] = arrayROC(x1, x2, reps, alpha, 'notpfp')
%
% DESCRIPTION
%   Computes the area under the Receiver Operating Characteristic (ROC)
%   curve for each corresponding column of two input matrices. The ROC
%   area (AUC) quantifies the discriminability between two distributions,
%   equivalent to the probability that a randomly drawn sample from x2
%   exceeds a randomly drawn sample from x1: AUC = P(X2 > X1).
%
%   An AUC of 0.5 indicates chance-level discrimination (overlapping
%   distributions), while 1.0 indicates perfect separation (all x2 > x1)
%   and 0.0 indicates perfect reversed separation (all x1 > x2).
%
% INPUTS
%   x1    - [nTrials1 x nBins] matrix of values for condition 1
%   x2    - [nTrials2 x nBins] matrix of values for condition 2
%           NaN values are automatically excluded on a per-bin basis.
%
% OPTIONAL INPUTS
%   reps  - Number of bootstrap replicates for confidence intervals.
%           Default: 200. Pass [] to skip CI computation.
%   alpha - Significance level for confidence intervals (default: 0.05)
%           Yields (1-alpha)*100% CIs, e.g., alpha=0.05 gives 95% CIs.
%   'notpfp' - Flag to suppress tpfp computation for faster execution.
%              When specified, tpfp output will be empty cell array.
%
% OUTPUTS
%   a     - [1 x nBins] vector of AUC values
%   ci    - [2 x nBins] matrix of confidence intervals [lower; upper]
%           Returns NaN if reps is empty or not specified.
%   tpfp  - {1 x nBins} cell array of ROC curve coordinates.
%           Each cell contains [nPoints x 2] matrix of [FPR, TPR] pairs.
%           Empty if 'notpfp' flag is set.
%   sig   - [1 x nBins] significance indicator vector:
%            +1 : AUC significantly above 0.5 (CI lower bound > 0.5)
%            -1 : AUC significantly below 0.5 (CI upper bound < 0.5)
%             0 : Not significant (CI contains 0.5) or CI not computed
%
% ALGORITHM
%   ROC curves are constructed using the histogram method:
%   1. Pool and sort values from both conditions
%   2. Define decision criteria at midpoints between adjacent values
%   3. Compute true positive rate (TPR) and false positive rate (FPR)
%      at each criterion via cumulative histograms
%   4. Integrate ROC curve using trapezoidal rule
%
%   Bootstrap CIs resample with replacement from each condition
%   independently, recomputing AUC for each replicate. This is
%   vectorized across all replicates for efficiency.
%
% EXAMPLES
%   % Basic usage - compute AUC only
%   auc = arrayROC(condition1_data, condition2_data);
%
%   % With bootstrap confidence intervals (500 replicates)
%   [auc, ci, ~, sig] = arrayROC(x1, x2, 500);
%
%   % Custom alpha level (99% CIs)
%   [auc, ci] = arrayROC(x1, x2, 1000, 0.01);
%
%   % Fast mode - skip ROC curve point computation
%   [auc, ci, ~, sig] = arrayROC(x1, x2, 500, 0.05, 'notpfp');
%
% NOTES
%   - Requires Parallel Computing Toolbox for parfor acceleration.
%     Falls back to serial execution if unavailable.
%   - Bins with fewer than 2 valid samples in either condition return
%     NaN for AUC and CI, and 0 for significance.
%   - Tie handling: ties are resolved by appropriate histogram binning,
%     yielding the standard AUC definition with ties.
%
% See also: perfcurve, ranksum, tiedrank

% -------------------------------------------------------------------------
% Input parsing and validation
% -------------------------------------------------------------------------
[reps, ciFlag, alphaVal, computeTpfp] = parseArgs(varargin{:});

% Validate input dimensions
if size(x1, 2) ~= size(x2, 2)
    error('arrayROC:DimensionMismatch', ...
        'x1 and x2 must have the same number of columns.');
end

nBins = size(x1, 2);

% -------------------------------------------------------------------------
% Handle NaN values - compute valid sample masks and counts per bin
% -------------------------------------------------------------------------
g1 = ~isnan(x1);  % [nTrials1 x nBins] logical mask
g2 = ~isnan(x2);  % [nTrials2 x nBins] logical mask
n1 = sum(g1, 1);  % [1 x nBins] valid count per bin for condition 1
n2 = sum(g2, 1);  % [1 x nBins] valid count per bin for condition 2

% -------------------------------------------------------------------------
% Preallocate outputs
% -------------------------------------------------------------------------
a   = zeros(1, nBins);
ci  = NaN(2, nBins);
sig = zeros(1, nBins);

if computeTpfp
    tpfp = cell(1, nBins);
else
    tpfp = {};
end

% -------------------------------------------------------------------------
% Process each bin in parallel
% -------------------------------------------------------------------------
parfor i = 1:nBins
    % Extract valid data for this bin
    x1_bin = x1(g1(:,i), i);
    x2_bin = x2(g2(:,i), i);
    
    % Compute ROC statistics for this bin
    [a(i), ci(:,i), tpfp_i, sig(i)] = processBin(...
        x1_bin, x2_bin, n1(i), n2(i), ...
        ciFlag, reps, alphaVal, computeTpfp);
    
    % Store tpfp if requested (parfor requires this structure)
    if computeTpfp
        tpfp{i} = tpfp_i;
    end
end

end


%% ========================================================================
%  PROCESSS SINGLE BIN
%  ========================================================================
function [a, ci, tpfp, sig] = processBin(x1, x2, n1, n2, ...
    ciFlag, reps, alphaVal, computeTpfp)
% Process a single column/bin of data
%
% Inputs:
%   x1, x2      - Column vectors of valid (non-NaN) values
%   n1, n2      - Sample sizes (length of x1, x2)
%   ciFlag      - Boolean: compute bootstrap CI?
%   reps        - Number of bootstrap replicates
%   alphaVal    - Significance level for CI
%   computeTpfp - Boolean: compute ROC curve points?
%
% Outputs:
%   a    - Scalar AUC value
%   ci   - [2x1] confidence interval [lower; upper]
%   tpfp - [nPoints x 2] ROC curve coordinates, or []
%   sig  - Significance indicator (-1, 0, or +1)

% --- Handle insufficient data ---
if n1 < 2 || n2 < 2
    a    = NaN;
    ci   = [NaN; NaN];
    tpfp = NaN;
    sig  = 0;
    return
end

% --- Handle empty or degenerate data ---
if isempty(x1) || isempty(x2)
    a    = NaN;
    ci   = [NaN; NaN];
    tpfp = NaN;
    sig  = 0;
    return
end

% -------------------------------------------------------------------------
% Build histogram bin edges (criteria) from pooled sorted data
% -------------------------------------------------------------------------
pooled = sort([x1; x2]);
diffs  = diff(pooled);

% Find minimum positive difference to set bin spacing
posDiffs = diffs(diffs > 0);

if isempty(posDiffs)
    % All values are identical - distributions perfectly overlap
    a    = 0.5;
    ci   = [0.5; 0.5];
    tpfp = [0 0; 1 1];
    sig  = 0;
    return
end

% Set bin edges at midpoints between adjacent unique values
% This ensures each unique value falls into exactly one bin
spacer = min(posDiffs) / 2;

% Build criteria vector for histc
% - Start below minimum value
% - Place edge after each data point
% - End above maximum value
if any(diffs == 0)
    % Ties present: use unique to avoid redundant bin edges
    linCrits = unique([pooled(1) - spacer; ...
                       pooled + spacer; ...
                       pooled(end) + 1]);
else
    % No ties: direct construction is faster
    linCrits = [pooled(1) - spacer; ...
                pooled + spacer; ...
                pooled(end) + 1];
end

% -------------------------------------------------------------------------
% Compute ROC curve via cumulative histograms
% -------------------------------------------------------------------------
% histc counts how many values fall into each bin
% Flip and cumsum to get cumulative rates from high to low threshold
h1 = histc(x1, linCrits);
h2 = histc(x2, linCrits);

% FPR = fraction of x1 (negatives) above threshold
% TPR = fraction of x2 (positives) above threshold
fpr = cumsum(flipud(h1)) / n1;
tpr = cumsum(flipud(h2)) / n2;

% Store ROC curve points if requested
if computeTpfp
    tpfp = [fpr, tpr];
else
    tpfp = [];
end

% Compute AUC via trapezoidal integration
a = computeAUC(fpr, tpr);

% -------------------------------------------------------------------------
% Bootstrap confidence intervals
% -------------------------------------------------------------------------
if ciFlag && reps > 0
    ci = bootstrapCI(x1, x2, n1, n2, linCrits, reps, alphaVal);
    
    % Determine significance based on CI relative to chance (0.5)
    if ci(1) > 0.5
        sig = 1;   % Significantly above chance
    elseif ci(2) < 0.5
        sig = -1;  % Significantly below chance
    else
        sig = 0;   % CI contains 0.5 (not significant)
    end
else
    ci  = [NaN; NaN];
    sig = 0;
end

end


%% ========================================================================
%  BOOTSTRAP CONFIDENCE INTERVAL
%  ========================================================================
function ci = bootstrapCI(x1, x2, n1, n2, linCrits, reps, alphaVal)
% Compute bootstrap confidence interval for AUC
%
% Uses the same histogram bin edges (linCrits) as the original AUC
% computation, enabling vectorized histc across all bootstrap replicates.
%
% Inputs:
%   x1, x2    - Original data vectors
%   n1, n2    - Sample sizes
%   linCrits  - Histogram bin edges from original computation
%   reps      - Number of bootstrap replicates
%   alphaVal  - Significance level (e.g., 0.05 for 95% CI)
%
% Output:
%   ci - [2x1] vector [lower_bound; upper_bound]

% Pre-generate all bootstrap sample indices
% Each column is one bootstrap replicate
idx1 = randi(n1, n1, reps);  % [n1 x reps]
idx2 = randi(n2, n2, reps);  % [n2 x reps]

% Resample data: result is [n x reps] matrices
x1_boot = x1(idx1);  % [n1 x reps]
x2_boot = x2(idx2);  % [n2 x reps]

% Vectorized histogram computation across all replicates
% histc operates column-wise when given a matrix
H1 = histc(x1_boot, linCrits);  % [nCrits x reps]
H2 = histc(x2_boot, linCrits);  % [nCrits x reps]

% Compute cumulative rates for all replicates simultaneously
% flipud + cumsum gives cumulative from high to low threshold
FPR = cumsum(flipud(H1), 1) / n1;  % [nCrits x reps]
TPR = cumsum(flipud(H2), 1) / n2;  % [nCrits x reps]

% Vectorized trapezoidal AUC across all bootstrap samples
% For each replicate: sum of (average height) * (width)
auc_boot = sum(((TPR(1:end-1,:) + TPR(2:end,:)) / 2) .* ...
               diff(FPR, 1, 1), 1);  % [1 x reps]

% Extract percentiles for confidence interval
ci = prctile(auc_boot, 100 * [alphaVal/2, 1 - alphaVal/2])';

end


%% ========================================================================
%  TRAPEZOIDAL AUC COMPUTATION
%  ========================================================================
function a = computeAUC(fpr, tpr)
% Compute area under ROC curve using trapezoidal rule
%
% Inputs:
%   fpr - False positive rate vector (x-coordinates)
%   tpr - True positive rate vector (y-coordinates)
%
% Output:
%   a - Scalar area under the curve
%
% The trapezoidal rule approximates the integral as sum of trapezoids:
%   Area = sum of (average height) * (width)
%        = sum of ((y[i] + y[i+1])/2) * (x[i+1] - x[i])

a = sum(((tpr(1:end-1) + tpr(2:end)) / 2) .* diff(fpr));

end


%% ========================================================================
%  INPUT ARGUMENT PARSING
%  ========================================================================
function [reps, ciFlag, alphaVal, computeTpfp] = parseArgs(varargin)
% Parse optional input arguments
%
% Outputs:
%   reps        - Number of bootstrap replicates (default: 200)
%   ciFlag      - Whether to compute CI (true if reps specified)
%   alphaVal    - Significance level (default: 0.05)
%   computeTpfp - Whether to compute ROC curve points (default: true)

% Defaults
reps        = 200;
alphaVal    = 0.05;
ciFlag      = false;
computeTpfp = true;

% Check for 'notpfp' flag anywhere in arguments
if any(strcmpi(varargin, 'notpfp'))
    computeTpfp = false;
    % Remove flag from varargin for subsequent processing
    varargin(strcmpi(varargin, 'notpfp')) = [];
end

% Parse remaining positional arguments
nArgs = length(varargin);

if nArgs >= 1
    if ~isempty(varargin{1})
        reps   = varargin{1};
        ciFlag = true;
    end
end

if nArgs >= 2
    alphaVal = varargin{2};
end

end