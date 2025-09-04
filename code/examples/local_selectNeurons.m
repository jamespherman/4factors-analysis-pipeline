function sc = local_selectNeurons(sc, gSac, giveFeed)
    giveFeed('local_selectNeurons: Identifying task-modulated neurons...');
    % ensure scSide is defined; default if not (though this should be set by local_defineScSide)
    if ~isfield(sc, 'scSide')
        giveFeed('WARNING in local_selectNeurons: sc.scSide not defined. Assuming "right" SC for T1=Contra logic.');
        sc.scSide = 'right'; % default assumption or could be an error condition
    end
    % determine contralateral target trials based on scSide
    if strcmpi(sc.scSide, 'right')      % Right SC: T1 (Left Target) is Contralateral
        selectedTrials = gSac.isHighSaliencyTrial & gSac.isT1Trial;
        giveFeed('local_selectNeurons: Right SC detected. Using T1 (Left Target) as Contralateral for neuron selection.');
    elseif strcmpi(sc.scSide, 'left') % Left SC: T2 (Right Target) is Contralateral
        selectedTrials = gSac.isHighSaliencyTrial & gSac.isT2Trial;
        giveFeed('local_selectNeurons: Left SC detected. Using T2 (Right Target) as Contralateral for neuron selection.');
    else % unknown SC side
        giveFeed('WARNING in local_selectNeurons: SC side is "unknown". Using T1 as default for neuron selection. This may not be optimal.');
        selectedTrials = gSac.isHighSaliencyTrial & gSac.isT1Trial; % default to T1 if side is unknown
    end
    trialInd = find(selectedTrials);
    if isempty(trialInd)
        giveFeed('WARNING in local_selectNeurons: No trials match selection criteria (e.g., HighSal & Contralateral). Skipping neuron selection.');
        % initialize fields to reflect no neurons selected / no data
        if sc.nClusters > 0
          sc.memSacCounts = zeros(sc.nClusters, 4, 0); % [nNeurons x nEpochs x nTrials]
          sc.sigEpochComparison = false(sc.nClusters, 3); % [nNeurons x nComparisons]
          sc.analysisNeurons = false(sc.nClusters, 1);    % [nNeurons x 1] logical
        else % no clusters to begin with
          sc.memSacCounts = [];
          sc.sigEpochComparison = [];
          sc.analysisNeurons = [];
        end
        return;
    end
    nTrials = numel(trialInd);
    nNeurons = sc.nClusters;
    if nNeurons == 0
        giveFeed('WARNING in local_selectNeurons: No neurons (sc.nClusters is 0). Skipping neuron selection.');
        sc.memSacCounts = zeros(0, 4, nTrials); % ensure correct empty dimensions
        sc.sigEpochComparison = false(0, 3);
        sc.analysisNeurons = false(0, 1);
        return;
    end
    memSacCounts = zeros(nNeurons, 4, nTrials); % to store firing rates for [Baseline, Visual, Delay, Saccade] epochs
    % define analysis windows relative to key trial events
    % windowsRel: {eventTimeVariable, offsetStart, offsetEnd, duration}
    % epochs: 1=Baseline, 2=Visual, 3=Delay, 4=Saccade
    windowsRelDefinition = { % This structure defines how windows are calculated for each trial
        'targetOnTime',  -0.075, 0.025,  0.1;   % Baseline: targetOn -0.075 to targetOn +0.025
        'targetOnTime',  0.05,   0.2,    0.15;  % Visual:   targetOn +0.05  to targetOn +0.2
        'fixOffTime',    -0.15,  0.05,   0.2;   % Delay:    fixOff   -0.15  to fixOff   +0.05
        'saccadeOnTime', -0.025, 0.05,   0.075  % Saccade:  saccadeOn-0.025 to saccadeOn+0.05
    };
    % iterate through selected trials to compute firing rates for different epochs
    for k_trial = 1:nTrials
        t = trialInd(k_trial); % actual trial index in gSac structure
        % iterate through defined epochs (Baseline, Visual, Delay, Saccade)
        for w_win = 1:size(windowsRelDefinition,1)
            eventNameStr = windowsRelDefinition{w_win,1};
            eventTimeVal = gSac.(eventNameStr)(t); % get event time for current trial
            
            startOffset  = windowsRelDefinition{w_win,2};
            endOffset    = windowsRelDefinition{w_win,3};
            winDur       = windowsRelDefinition{w_win,4};
            if isnan(eventTimeVal)
                giveFeed(['Warning: NaN event time for ' eventNameStr ' on trial index ' num2str(t) '. Setting epoch ' num2str(w_win) ' FR to NaN.']);
                memSacCounts(:, w_win, k_trial) = NaN;
                continue; % to next epoch for this trial if event time is NaN
            end
            
            startTime = eventTimeVal + startOffset;
            endTime   = eventTimeVal + endOffset;
            
            memSacCounts(:, w_win, k_trial) = local_computeFiringRate(sc, startTime, endTime, winDur, nNeurons);
        end
    end
    sc.memSacCounts = memSacCounts; % store epoch-based firing rates
    % perform statistical tests to identify task-modulated neurons
    sigEpochComparison = false(nNeurons, 3); % [Visual vs Base, Delay vs Base, Saccade vs Base]
    analysisNeurons = false(nNeurons, 1);    % logical vector indicating selected neurons
    if nTrials < 2 % not enough trials for statistical comparison
        giveFeed('Warning: Too few trials for statistical neuron selection. Marking all neurons as not analyzed.');
        sc.sigEpochComparison = sigEpochComparison;
        sc.analysisNeurons = analysisNeurons;
        return;
    end
    % iterate through neurons to test for significant modulation across epochs
    for i_neuron = 1:nNeurons
        X = squeeze(memSacCounts(i_neuron, :, :))'; % X is [nTrials x nEpochs] for the current neuron
        X = X(~any(isnan(X),2),:); % remove trials with NaN in any epoch for this neuron
    
        if size(X,1) < 2 % not enough valid trials for this neuron
            continue;
        end
        
        try
            % use a non-parametric test (Friedman) to check for overall epoch modulation
            % X is [nTrials x nEpochs]. The second argument '1' specifies one
            % observation per group. 'off' suppresses the table display.
            [p_friedman, ~] = friedman(X, 1, 'off'); 

            % if the overall test is significant, proceed to pairwise comparisons
            if p_friedman < 0.05
                alpha_corr = 0.05 / 3; % bonferroni correction for 3 comparisons
                colMap = [1 2; 1 3; 1 4]; % compare Visual(2), Delay(3), Saccade(4) against Baseline(1)
                
                for j_comp = 1:size(colMap, 1)
                    % use ROC AUC to test if epoch activity is significantly different from baseline.
                    % NOTE: This assumes 'arrayROC' is a valid, existing helper function.
                    [~, ci] = arrayROC(X(:, colMap(j_comp, 1)), X(:, colMap(j_comp, 2)), 200, alpha_corr);
                    
                    % check if the 95% confidence interval for the AUC excludes 0.5
                    sigEpochComparison(i_neuron, j_comp) = (ci(1) > 0.5 || ci(2) < 0.5);
                end
            end
        catch ME_stats_select
            giveFeed(['Statistical test (Friedman or ROC) failed for neuron ' num2str(i_neuron) ': ' ME_stats_select.message]);
        end
        
        % select neuron if it shows significant modulation in any epoch and has a decent firing rate
        if any(sigEpochComparison(i_neuron, :)) && (nTrials==0 || max(mean(X,1,'omitnan')) > 5 ) % check mean FR across epochs
            analysisNeurons(i_neuron) = true;
        end
    end
    % broaden selection criteria if too few neurons are initially selected
    if nnz(analysisNeurons) < min(10, nNeurons) && nNeurons > 0
        giveFeed(['Too few neurons initially selected (' num2str(nnz(analysisNeurons)) '). Broadening criteria to any sig epoch comparison...']);
        analysisNeurons = any(sigEpochComparison,2); % include if any epoch comparison was significant
    end
    
    if nNeurons == 0 % handle case of no neurons
        analysisNeurons = false(0,1);
    end
    sc.sigEpochComparison = sigEpochComparison;
    sc.analysisNeurons = analysisNeurons; % store the final list of selected neurons
end