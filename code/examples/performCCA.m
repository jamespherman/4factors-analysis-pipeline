function ccaResults = performCCA(scMatrices, sncMatrices, varargin)
% performs cca across sc and snc simultaneous recordings

% parse opts
p = inputParser;
addParameter(p, 'taskType', 'gSac', @ischar);
parse(p, varargin{:});
taskType = p.Results.taskType;

% set task-specific params
if strcmpi(taskType, 'gSac')
    binWidth = 0.02; % s
    preTimes = [-0.1, -0.5, -0.1];
    postTimes = [0.9, 0.5, 0.9];
    timeMap.targetOn   = 1;
    timeMap.saccadeOn  = 2;
    timeMap.reward     = 3;
    analysisWins.targetOn   = [0.05,  0.55];
    analysisWins.saccadeOn  = [-0.25,  0.25];
    analysisWins.reward     = [0.1,  0.6];
    conds = {'lowSal_lowRew_t1','lowSal_lowRew_t2','lowSal_highRew_t1','lowSal_highRew_t2',...
        'highSal_lowRew_t1','highSal_lowRew_t2','highSal_highRew_t1','highSal_highRew_t2'};
elseif strcmpi(taskType, 'attn')
    binWidth = 0.02; % s
    preTimes = [-0.1, -0.1, -0.1, -0.1];
    postTimes = [0.9, 0.9, 0.9, 0.9];
    timeMap.cueOn   = 1;
    timeMap.stimOn  = 2;
    timeMap.stimChg = 3;
    timeMap.reward  = 4;
    analysisWins.cueOn   = [0.05, 0.55];
    analysisWins.stimOn  = [0.05, 0.55];
    analysisWins.stimChg = [-0.25, 0.25];
    analysisWins.reward  = [0.1, 0.6];

    conds = {'cueIn_chgIn_hit','cueIn_chgIn_miss','cueOut_chgOut_hit','cueOut_chgOut_miss'};
else
    error('unknown taskType. use ''gSac'' or ''attn''.');
end

%% params
% loop over time windows
timeFields = fieldnames(analysisWins);
for iTF = 1:numel(timeFields)
    tw = timeFields{iTF};

    % build time vector for the full window (50 bins)
    idx = timeMap.(tw);
    tVec = linspace(preTimes(idx), postTimes(idx), 50);

    % find bin indices for the sub-window we want
    subWin = analysisWins.(tw);
    binMask = (tVec >= subWin(1)) & (tVec <= subWin(2));

    %% now do each condition
    for c = 1:numel(conds)
        condName = conds{c};
        
        % For the reward window in attn, skip miss conditions because
        % misses do not yield a reward. Just assign 0.
        if strcmpi(taskType, 'attn') && strcmp(tw, 'reward') && contains(condName, 'miss')
            ccaResults.(tw).(condName).R = 0;
            continue;
        end

        % sc data
        scData = scMatrices.(tw).(condName); % [nNeurons x nTrials x nBins]
        % snc data
        sncData = sncMatrices.(tw).(condName);

        % average across bins in sub-window (if needed)
        scAvg = mean(scData(:,:,binMask), 3); % [nNeurons x nTrials]
        sncAvg = mean(sncData(:,:,binMask), 3);

        % instead of averaging, rearrange time-points to yield (nBins X nTrials) X nNeurons
        scRerng = reshape(permute(scData(:,:,binMask), [3 2 1]), ...
            size(scData, 2) * nnz(binMask), size(scData, 1));
        sncRerng = reshape(permute(sncData(:,:,binMask), [3 2 1]), ...
            size(sncData, 2) * nnz(binMask), size(sncData, 1));

        % Check if there are at least 2 observations and at least one non-constant column in both matrices
        if size(scRerng, 1) < 2 || size(sncRerng, 1) < 2 || ...
                all(var(scRerng,0,1)==0) || all(var(sncRerng,0,1)==0)
            fprintf('Skipping condition %s in time window %s due to insufficient data.\n', condName, tw);
            continue; % Skip this condition if not enough data
        end

        % Perform canonical correlation analysis
        [A, B, R, U, V, stats] = canoncorr(scRerng, sncRerng);
        % R is vector of canonical correlations
        % stats contains significance tests for each dimension

        % Store Correlations (Loadings) for All Canonical Dimensions
        for d = 1:length(R)
            Ud = U(:,d);
            Vd = V(:,d);
            corrXUd = corr(scRerng, Ud);
            corrYVd = corr(sncRerng, Vd);

            % Store them in a cell or matrix
            ccaResults.(tw).(condName).varPortionSC(:, d) = corrXUd.^2;
            ccaResults.(tw).(condName).varPortionSNc(:, d) = corrYVd.^2;
        end

        % store immediate results
        ccaResults.(tw).(condName).A = A;
        ccaResults.(tw).(condName).B = B;
        ccaResults.(tw).(condName).R = R;
        ccaResults.(tw).(condName).U = U;
        ccaResults.(tw).(condName).V = V;
        ccaResults.(tw).(condName).stats = stats;
    end
end