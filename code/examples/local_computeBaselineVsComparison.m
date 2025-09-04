function baselineComp = local_computeBaselineVsComparison(rawArray, binCenters, preTime, postTime, eventName, condContra, condIpsi, giveFeed)
    % computes ROC AUC comparing baseline activity to post-event activity for contra and ipsi conditions.
    % rawArray: [nTrials x nNeurons x nBins]
    % binCenters: [1 x nBins]
    % condContra, condIpsi: logical vectors for trial selection
    giveFeed('local_computeBaselineVsComparison: Computing baseline vs. comparison discriminability...');

    % define baseline and comparison (post-event) windows based on binCenters
    if strcmp(eventName, 'sac_saccadeOn') % special baseline for saccade-aligned activity
         baselineWindow = [-0.5, -0.4]; % specific pre-saccadic window
    else % general baseline window (typically pre-event)
         baselineWindow = [preTime, 0]; % assumes preTime is negative (e.g., -0.25 to 0)
         if preTime >=0 % handle cases where preTime might not be suitable for baseline
             baselineWindow = [-0.2, 0]; % default pre-event baseline
             giveFeed(['Warning: preTime >=0 for ' eventName ', using default baseline [-0.2, 0] for baselineVsComp.']);
         end
    end
    compWindow = [0, postTime]; % comparison window is from event onset to postTime
    if postTime <=0 % handle cases where postTime might not be suitable
        compWindow = [0, 0.5]; % default post-event window
        giveFeed(['Warning: postTime <=0 for ' eventName ', using default comp window [0, 0.5] for baselineVsComp.']);
    end

    % get indices for baseline and comparison bins
    baselineIdx = (binCenters >= baselineWindow(1)) & (binCenters < baselineWindow(2));
    compIdx     = (binCenters >= compWindow(1)) & (binCenters <= compWindow(2));

    % check if any bins fall within the defined windows
    if ~any(baselineIdx) || ~any(compIdx)
        giveFeed(['Warning in local_computeBaselineVsComparison for ' eventName ': No bins in baseline or comparison window. BaselineWin: [' num2str(baselineWindow) '], CompWin: [' num2str(compWindow) ']. Returning empty.']);
        baselineComp = struct('contra', struct('roc',[],'sig',[]), 'ipsi', struct('roc',[],'sig',[]), 'combined', struct('roc',[],'sig',[]));
        return;
    end

    baselineComp = struct(); % initialize results structure

    % nested function to compute ROC for a given trial mask (e.g., condContra)
    function [rocVals, sigVals] = computeROCForMaskLocal(trialMask)
        % trialMask is a logical vector for selecting trials from rawArray's first dimension
        if sum(trialMask) < 2 % need at least 2 trials for ROC
            rocVals = NaN(size(rawArray,2), sum(compIdx)); % [nNeurons x nComparisonBins]
            sigVals = NaN(size(rawArray,2), sum(compIdx));
            return;
        end
        
        subData = rawArray(trialMask, :, :); % select trials: [nSelectedTrials x nNeurons x nAllBins]
        [nSelTrials, nNeuronsSub, ~] = size(subData);
        rocVals = NaN(nNeuronsSub, sum(compIdx)); % initialize results for this mask
        sigVals = NaN(nNeuronsSub, sum(compIdx));
        
        % iterate through each neuron
        for iNeuron = 1:nNeuronsSub
            % get baseline activity for this neuron: average FR across baseline bins, for each selected trial
            neuronBaselineDataPerTrial = mean(squeeze(subData(:, iNeuron, baselineIdx)), 2, 'omitnan'); % result: [nSelectedTrials x 1]
            
            % get comparison activity: FR in each comparison bin, for each selected trial
            neuronCompDataPerTrialPerBin = squeeze(subData(:, iNeuron, compIdx)); % result: [nSelectedTrials x nComparisonBins]
            
            % handle cases of single trial or single comparison bin after squeeze
            if nSelTrials == 1 && sum(compIdx) > 1 % if 1 trial, multiple comp bins
                neuronCompDataPerTrialPerBin = neuronCompDataPerTrialPerBin'; % ensure [1 x nComparisonBins]
            end
            if isempty(neuronBaselineDataPerTrial) || isempty(neuronCompDataPerTrialPerBin) || nSelTrials < 2
                 rocVals(iNeuron,:) = NaN; sigVals(iNeuron,:) = NaN; % not enough data
                continue;
            end
            
            currentNeuronRoc = NaN(1, sum(compIdx));
            currentNeuronSig = NaN(1, sum(compIdx));
            % compare baseline distribution to distribution in each comparison bin
            for iBinComp = 1:sum(compIdx)
                compDataThisBin = neuronCompDataPerTrialPerBin(:, iBinComp); % [nSelectedTrials x 1] for current comparison bin
                % ensure enough non-NaN data points for ROC
                if sum(~isnan(neuronBaselineDataPerTrial)) <2 || sum(~isnan(compDataThisBin)) < 2
                    continue;
                end
                try
                    % arrayROC compares two distributions (baseline vs current comparison bin)
                    [roc_val, ~, ~, sig_val] = arrayROC(neuronBaselineDataPerTrial, compDataThisBin, 200);
                    currentNeuronRoc(iBinComp) = roc_val;
                    currentNeuronSig(iBinComp) = sig_val;
                catch ME_roc_base
                    giveFeed(['arrayROC failed in baselineVsComp (' eventName ') for neuron ' num2str(iNeuron) ', bin ' num2str(iBinComp) ': ' ME_roc_base.message]);
                end
            end
            rocVals(iNeuron,:) = currentNeuronRoc;
            sigVals(iNeuron,:) = currentNeuronSig;
        end
    end % end of computeROCForMaskLocal nested function

    % compute for contralateral trials
    [rocContra, sigContra] = computeROCForMaskLocal(condContra);
    baselineComp.contra.roc = rocContra; baselineComp.contra.sig = sigContra;

    % compute for ipsilateral trials
    [rocIpsi, sigIpsi] = computeROCForMaskLocal(condIpsi);
    baselineComp.ipsi.roc = rocIpsi; baselineComp.ipsi.sig = sigIpsi;

    % compute for combined (contra or ipsi) trials
    [rocCombined, sigCombined] = computeROCForMaskLocal(condContra | condIpsi);
    baselineComp.combined.roc = rocCombined; baselineComp.combined.sig = sigCombined;
end