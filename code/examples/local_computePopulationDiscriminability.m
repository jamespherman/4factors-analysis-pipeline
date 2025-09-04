function popDiscrim = local_computePopulationDiscriminability(dataArray, cond1, cond2, giveFeed)
    % computes population discriminability using SVM for each time bin.
    % dataArray: [nTrials x nNeurons x nTimeBins]
    % cond1, cond2: logical trial selectors
    giveFeed('local_computePopulationDiscriminability: Computing population SVM discriminability...');
    popDiscrim = struct(); % initialize results
    [~, nNeurons, nTimeBins] = size(dataArray);
    
    % initialize population discriminability fields
    popDiscrim.accuracy = NaN(1, nTimeBins); % [1 x nTimeBins]
    popDiscrim.ci       = NaN(2, nTimeBins); % [2 x nTimeBins] for CI bounds

    minTrialsForSVM = 5; % minimum trials per condition needed for SVM

    if nNeurons == 0
        giveFeed('Warning: No neurons for population discriminability.');
        return; % fields remain NaN
    end

    % iterate through each time bin
    for iBin = 1:nTimeBins
        % extract data for the current bin: [nTrials x nNeurons]
        binData1_allNeurons = squeeze(dataArray(cond1, :, iBin));
        binData2_allNeurons = squeeze(dataArray(cond2, :, iBin));
        
        % handle squeeze behavior for single neuron or single trial cases
        if nNeurons == 1 && sum(cond1) > 1, binData1_allNeurons = binData1_allNeurons(:); end % ensure column vector
        if nNeurons == 1 && sum(cond2) > 1, binData2_allNeurons = binData2_allNeurons(:); end % ensure column vector
        if sum(cond1) == 1 && nNeurons > 1, binData1_allNeurons = binData1_allNeurons(:)'; end % ensure row vector
        if sum(cond2) == 1 && nNeurons > 1, binData2_allNeurons = binData2_allNeurons(:)'; end % ensure row vector

        numCond1Trials = size(binData1_allNeurons, 1);
        numCond2Trials = size(binData2_allNeurons, 1);

        % check if enough trials for SVM
        if numCond1Trials < minTrialsForSVM || numCond2Trials < minTrialsForSVM
            % popDiscrim fields already NaN for this bin
            continue;
        end
        
        % prepare data (X) and labels (Y) for SVM
        X = [binData1_allNeurons; binData2_allNeurons]; % [totalSelectedTrials x nNeurons]
        Y = [zeros(numCond1Trials, 1); ones(numCond2Trials, 1)]; % labels: 0 for cond1, 1 for cond2
        
        % remove rows with NaNs (e.g., if a neuron had NaN FR for a trial in this bin)
        nanRows = any(isnan(X),2);
        X(nanRows,:) = [];
        Y(nanRows,:) = [];
        
        % check again after NaN removal if enough samples remain
        if size(X,1) < (2*minTrialsForSVM) || isempty(X)
            continue;
        end

        try
            % train SVM classifier with cross-validation
            svmObj = fitcsvm(X, Y, 'Standardize', true, 'KernelFunction', 'linear', 'Prior', 'uniform', 'CrossVal', 'on', 'KFold', 5);
            yHat = kfoldPredict(svmObj); % get cross-validated predictions
            
            correctPredictions = sum(yHat == Y);
            totalPredictions = length(Y);
            if totalPredictions > 0
                [p, pCI] = binofit(correctPredictions, totalPredictions); % binomial fit for accuracy and CI
                popDiscrim.accuracy(iBin) = p;
                popDiscrim.ci(:, iBin)    = pCI';
            end
        catch ME_svm
            giveFeed(['SVM failed at bin ' num2str(iBin) ': ' ME_svm.message]);
            % values for this bin remain NaN
        end
    end
end