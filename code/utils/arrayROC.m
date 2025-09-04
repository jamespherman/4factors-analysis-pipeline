function [a, ci, tpfp, sig] = arrayROC(x1, x2, varargin)
% [a, ci, tpfp, sig] = arrayROC(x1, x2, [reps for bootstrapped ci], alpha)
%
% Returns:
%   a   - Area under ROC curve for each time bin
%   ci  - 95% confidence intervals [2 x nBins]
%   tpfp - TPR/FPR points for ROC curve
%   sig  - Significance vector: +1 (above 0.5), -1 (below 0.5), 0 
%           (not significant)

% Parse input arguments
reps = 200;


switch nargin
    case 3
        reps = varargin{1};
        ciFlag = true;
        alphaVal = 0.05;
    case 4
        if ~isempty(varargin{1})
            reps = varargin{1};
            ciFlag = true;
        else
            ciFlag = false;
        end
        alphaVal = varargin{2};
    otherwise
        ciFlag = false;
        reps = 0;
        alphaVal = 0.05;
end

% Get dimensions
[nTrials1, nBins] = size(x1);
nTrials2 = size(x2, 1);

% Preallocate outputs
a = zeros(1, nBins);
ci = zeros(2, nBins);
tpfp = cell(1, nBins);
sig = zeros(1, nBins);  % New significance output

% Handle NaNs
g1 = ~isnan(x1);
g2 = ~isnan(x2);
n1 = sum(g1);
n2 = sum(g2);

% Preallocate temporary arrays for parfor
temp_a = zeros(1, nBins);
temp_ci = zeros(2, nBins);
temp_tpfp = cell(1, nBins);
temp_sig = zeros(1, nBins);

try
    % Loop over columns
    parfor i = 1:nBins
        % Get valid data for this bin
        x1_bin = x1(g1(:,i), i);
        x2_bin = x2(g2(:,i), i);
        
        % Skip if insufficient samples
        if n1(i) < 2 || n2(i) < 2
            temp_a(i) = NaN;
            temp_tpfp{i} = NaN;
            temp_ci(:,i) = NaN;
            temp_sig(i) = 0;
            continue;
        end
        
        % Create input for vectorROC
        vrInput = sortrows([[x1_bin; x2_bin], ...
            [zeros(n1(i), 1); ones(n2(i), 1)]]);
        
        % Compute ROC statistics
        if ciFlag
            [temp_a(i), temp_tpfp{i}, this_ci] = vectorROC(vrInput, ...
                n1(i), n2(i), reps, alphaVal);
            temp_ci(:,i) = this_ci;
            
            % Determine significance
            if this_ci(1) > 0.5
                temp_sig(i) = 1;      % Significantly above 0.5
            elseif this_ci(2) < 0.5
                temp_sig(i) = -1;     % Significantly below 0.5
            end
        else
            [temp_a(i), temp_tpfp{i}] = vectorROC(vrInput, n1(i), n2(i));
            temp_ci(:,i) = NaN;
        end
    end
    
    % Copy temporary arrays to output
    a = temp_a;
    ci = temp_ci;
    tpfp = temp_tpfp;
    sig = temp_sig;
    
catch me
    keyboard
end

end

function [a, tpfp, varargout]  = vectorROC(Z, n1, n2, varargin)

% parse Z into x1/x2:
x1 = Z(Z(:,2) == 0, 1);
x2 = Z(Z(:,2) == 1, 1);

spacer      = mingz(diff(Z(:,1)))/2;

% Parse optional arguments
if nargin > 3 && ~isempty(varargin)
    reps = varargin{1};
end
if nargin > 4 && length(varargin) > 1
    alphaVal = varargin{2};
end

% if the data are weird, return something anyway.
if isempty(spacer)
    a = 0.5;
    tpfp = [0 0; 1 1];
    if nargout > 2
        varargout{1} = [0.5; 0.5];
    end
else
    if nnz(diff(Z(:,1))==0)>0
        linCrits = unique([[Z(1,1); Z(:,1)] + ...
            spacer*[-1; ones(n1+n2, 1)]; Z(end,1) + 1]);
        nfb = length(linCrits);
    else
        linCrits = [[Z(1,1); Z(:,1)] + ...
            spacer*[-1; ones(n1+n2, 1)]; Z(end,1) + 1];
        nfb = n1 + n2 + 2;
    end
    try
        tpfp        = [cumsum(flipud(histc(x1, linCrits)))/n1, ...
            cumsum(flipud(histc(x2, linCrits)))/n2];
    catch me
        keyboard
    end
    a           = auc(tpfp);

    % compute Bootstrapped CI if desired
    if nargout > 2
        TPFP(:,1,:) = cumsum(reshape(flipud(...
            histc(x1(randi(n1,n1,reps)), linCrits)), nfb,1,reps))/n1;
        TPFP(:,2,:) = cumsum(reshape(flipud(...
            histc(x2(randi(n2,n2,reps)), linCrits)), nfb,1,reps))/n2;
        as = auc(TPFP);
        varargout{1} = prctile(as, 100*[(alphaVal/2) (1 - alphaVal/2)]);
    end
end

if nargout > 2
    if all(isnan(varargout{1}(:)))
        keyboard
    end
end

end

function a = auc(XY)
a = sum(squeeze(mean([XY(1:end-1,2,:), XY(2:end,2,:)],2) .* ...
    diff(XY(:,1,:))));
end

function y = mingz(x)
y = min(x(x>0));
end