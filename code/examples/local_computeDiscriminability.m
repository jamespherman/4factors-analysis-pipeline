function discrim = local_computeDiscriminability(dataArray, cond1, cond2, giveFeed)
    % computes single neuron ROC AUC and population SVM discriminability between two conditions.
    % dataArray: [nTrials x nSelectedNeurons x nTimeBins]
    % cond1, cond2: logical vectors for trial selection under each condition
    giveFeed('local_computeDiscriminability: Computing single neuron and population discriminability...');
    discrim = struct(); % initialize results structure
    [~, nNeurons, nTimeBins] = size(dataArray);

    if nNeurons == 0
        giveFeed('Warning: No neurons in dataArray for discriminability. Returning empty.');
        discrim.singleNeuron.auc = []; discrim.singleNeuron.ci = []; discrim.singleNeuron.sig = [];
        discrim.population.accuracy = []; discrim.population.ci = [];
        return;
    end

    % initialize single neuron discriminability fields
    discrim.singleNeuron.auc = NaN(nNeurons, nTimeBins); % [nNeurons x nTimeBins]
    discrim.singleNeuron.ci  = NaN(2, nNeurons, nTimeBins); % [2_x_nNeurons_x_nTimeBins] for CI bounds
    discrim.singleNeuron.sig = NaN(nNeurons, nTimeBins); % significance (p-value or logical)

    % select data for the two conditions
    data1_trials = dataArray(cond1, :, :); % [nCond1Trials x nNeurons x nTimeBins]
    data2_trials = dataArray(cond2, :, :); % [nCond2Trials x nNeurons x nTimeBins]

    % check if enough trials per condition for ROC
    if size(data1_trials,1) < 2 || size(data2_trials,1) < 2
        giveFeed('Warning: Less than 2 trials for one or both conditions in local_computeDiscriminability. Single neuron AUCs will be NaN.');
        % fields are already NaN initialized
    else
        % iterate through each neuron
        for iNeuron = 1:nNeurons
            neuronData1 = squeeze(data1_trials(:, iNeuron, :)); % [nCond1Trials x nTimeBins]
            neuronData2 = squeeze(data2_trials(:, iNeuron, :)); % [nCond2Trials x nTimeBins]
            
            % handle squeeze behavior if only 1 trial or 1 time bin
            if size(data1_trials,1) == 1 && nTimeBins > 1, neuronData1 = neuronData1'; end % ensure [1 x nTimeBins]
            if size(data2_trials,1) == 1 && nTimeBins > 1, neuronData2 = neuronData2'; end % ensure [1 x nTimeBins]
            if nTimeBins == 1 && size(data1_trials,1) > 1, neuronData1 = neuronData1(:); end % ensure [nTrials x 1]
            if nTimeBins == 1 && size(data2_trials,1) > 1, neuronData2 = neuronData2(:); end % ensure [nTrials x 1]

            if isempty(neuronData1) || isempty(neuronData2) % skip if no data for this neuron
                continue;
            end
            
            try
                % arrayROC computes ROC AUC for each time bin by comparing neuronData1 vs neuronData2
                [auc, ci, ~, sig] = arrayROC(neuronData1, neuronData2, 200);
                discrim.singleNeuron.auc(iNeuron, :)   = auc; % [1 x nTimeBins]
                discrim.singleNeuron.ci(:, iNeuron, :) = ci;  % [2 x nTimeBins] -> store as [2 x 1 x nTimeBins] then reshape
                discrim.singleNeuron.sig(iNeuron, :)   = sig; % [1 x nTimeBins]
            catch ME_roc_discrim
                giveFeed(['arrayROC failed in local_computeDiscriminability for neuron ' num2str(iNeuron) ': ' ME_roc_discrim.message]);
                % values will remain NaN
            end
        end
    end
    
    % compute population discriminability using SVM
    discrim.population = local_computePopulationDiscriminability(dataArray, cond1, cond2, giveFeed);
end