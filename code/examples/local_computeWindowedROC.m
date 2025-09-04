function sc_out = local_computeWindowedROC(sc_in, eventName, giveFeed)
    % computes ROC AUC for defined comparisons within a specific time window for an event.
    % uses sc.taskEventTimeWindows to define the analysis window.
    giveFeed(['local_computeWindowedROC: Computing for event: ' eventName]);
    sc_out = sc_in; % operate on a copy, assign back if needed

    % check if a specific analysis time window is defined for this event
    if ~isfield(sc_out, 'taskEventTimeWindows') || ~isfield(sc_out.taskEventTimeWindows, eventName)
        giveFeed(['No time window defined in sc.taskEventTimeWindows for event: ' eventName '. Skipping windowed ROC.']);
        return;
    end
    cfg = sc_out.analysisConfig.(eventName); % config for the current event

    % ensure necessary data arrays are present
    if ~isfield(cfg, 'rawArray') || ~isfield(cfg, 'binCenters') || ~isfield(cfg, 'condVars')
        giveFeed(['rawArray, binCenters, or condVars missing for ' eventName '. Skipping windowed ROC.']);
        return;
    end

    % get the pre-defined time window [startTime, endTime] for this event
    timeWindow = sc_out.taskEventTimeWindows.(eventName);
    % find bin indices that fall within this time window
    binIndices = (cfg.binCenters >= timeWindow(1)) & (cfg.binCenters <= timeWindow(2));

    if sum(binIndices) == 0 % no bins in the specified window
        giveFeed(sprintf('No time bins in window [%g, %g] for %s. Skipping windowed ROC.', timeWindow(1), timeWindow(2), eventName));
        return;
    end

    % initialize results structure for windowed ROC
    cfg.windowedROC = struct();
    cfg.windowedROC.timeWindow = timeWindow;
    cfg.windowedROC.binIndices = binIndices; % store which bins were used

    % determine which set of comparisons to use (saccade or attention task)
    if startsWith(eventName, 'sac_')
        if isfield(sc_out.comparisons, 'gSac')
            comparisons = sc_out.comparisons.gSac;
        else
            giveFeed(['Warning: sc_out.comparisons.gSac not found for ' eventName '. Skipping windowed ROC.']); return;
        end
    else % assume 'attn_' prefix
        if isfield(sc_out.comparisons, 'attn')
            comparisons = sc_out.comparisons.attn;
        else
             giveFeed(['Warning: sc_out.comparisons.attn not found for ' eventName '. Skipping windowed ROC.']); return;
        end
    end

    % iterate through defined comparisons for this task type
    for iComp = 1:length(comparisons)
        comp = comparisons(iComp); % current comparison (e.g., targetLoc)
        rawCompName = ['raw_' comp.name]; % name for storing results, matches discriminability field names
        
        cond1Field = comp.condition1; cond2Field = comp.condition2;
        % ensure condition vectors exist in condVars
        if ~isfield(cfg.condVars, cond1Field) || ~isfield(cfg.condVars, cond2Field)
            giveFeed(['CondVars fields ' cond1Field ' or ' cond2Field ' missing for ' comp.name ' in ' eventName '. Skipping windowed ROC for this comparison.']); continue;
        end
        cond1 = cfg.condVars.(cond1Field); % logical vector for condition 1 trials
        cond2 = cfg.condVars.(cond2Field); % logical vector for condition 2 trials

        % check for sufficient trials
        if sum(cond1) < 2 || sum(cond2) < 2
            giveFeed(['Not enough trials for ' comp.name ' for windowed ROC in ' eventName '. Cond1 trials: ' num2str(sum(cond1)) ', Cond2 trials: ' num2str(sum(cond2)) '. Skipping.']); continue;
        end
        
        rawArray = cfg.rawArray; % [nAllTrialsForEvent x nSelectedNeurons x nAllBins]
        [~, nNeurons, ~] = size(rawArray);
        if nNeurons == 0, giveFeed(['No selected neurons for ' eventName '. Skipping windowed ROC for ' comp.name]); continue; end
        
        % calculate mean activity within the window for each trial and neuron
        meanActivity1_trials_neurons = NaN(sum(cond1), nNeurons); % [nCond1Trials x nNeurons]
        meanActivity2_trials_neurons = NaN(sum(cond2), nNeurons); % [nCond2Trials x nNeurons]
        
        % extract data for condition 1 trials within the window
        data_cond1_window = rawArray(cond1, :, binIndices); % [nCond1Trials x nNeurons x nBinsInWindow]
        % extract data for condition 2 trials within the window
        data_cond2_window = rawArray(cond2, :, binIndices); % [nCond2Trials x nNeurons x nBinsInWindow]

        % average across bins in the window for each neuron and trial
        for n_neuron_idx = 1:nNeurons % iterate over selected neurons
            meanActivity1_trials_neurons(:, n_neuron_idx) = mean(squeeze(data_cond1_window(:, n_neuron_idx, :)), 2, 'omitnan');
            meanActivity2_trials_neurons(:, n_neuron_idx) = mean(squeeze(data_cond2_window(:, n_neuron_idx, :)), 2, 'omitnan');
        end
        
        % initialize ROC results for each neuron
        aucValues = NaN(nNeurons, 1); ciValues = NaN(2, nNeurons); sigValues = NaN(nNeurons, 1);

        % compute ROC for each neuron using the window-averaged activity
        for n_neuron_idx = 1:nNeurons
            activity_N_C1 = meanActivity1_trials_neurons(:, n_neuron_idx); % [nCond1Trials x 1]
            activity_N_C2 = meanActivity2_trials_neurons(:, n_neuron_idx); % [nCond2Trials x 1]
            
            % remove NaNs that might result if all bins in window were NaN for a trial
            activity_N_C1 = activity_N_C1(~isnan(activity_N_C1));
            activity_N_C2 = activity_N_C2(~isnan(activity_N_C2));

            % ensure enough valid trial data for ROC computation
            if length(activity_N_C1) < 2 || length(activity_N_C2) < 2
                continue;
            end
            try
                [auc, ci, ~, sig] = arrayROC(activity_N_C1, activity_N_C2, 200, false); % 'false' for no figure
                aucValues(n_neuron_idx) = auc;
                ciValues(:, n_neuron_idx) = ci; % [2x1]
                sigValues(n_neuron_idx) = sig;
            catch ME_roc_win
                 giveFeed(['arrayROC failed in local_computeWindowedROC for ' eventName ', comp ' comp.name ', neuron ' num2str(n_neuron_idx) ': ' ME_roc_win.message]);
            end
        end
        
        % store results for this comparison
        rocResults = struct('auc', aucValues, 'ci', ciValues, 'sig', sigValues, ...
                            'meanAuc', mean(aucValues(~isnan(aucValues))), ... % mean AUC across neurons
                            'conditionLabels', {{comp.condition1, comp.condition2}});
        cfg.windowedROC.(rawCompName) = rocResults;
    end % end comparisons loop
    sc_out.analysisConfig.(eventName) = cfg; % update the config in the main structure
end