function data = local_defineAnalysisConfig(g4f, attn, data, giveFeed)
    giveFeed('local_defineAnalysisConfig: defining analysis event configurations...');
    comps = data.comparisons;          % pre-defined experimental comparisons (4f + attn)
    analysisConfig = struct();         % holds configurations per event

    %% gSac (4-factor) configuration
    % event list and windows (target-aligned for stim/prob/sal; saccade-aligned for movement/reward)
    gSacEventLabels = {'fixOn', 'targetOn', 'fixOff', 'saccadeOn', 'reward'};
    gSacPreTimes    = [-0.25,   -0.10,      -0.25,    -0.50,       -0.10];
    gSacPostTimes   = [ 0.75,    0.90,       0.75,     0.50,        0.90];

    % collect unique condition fields from 4f gSac comparisons
    if isempty(comps.gSac)
        gSacConditionFields = {};
        giveFeed('warning: no gSac comparisons defined. gSac analysis will be limited.');
    else
        nGsacComps = numel(comps.gSac);
        tmp = cell(1, nGsacComps * 2);
        for i_comp = 1:nGsacComps
            tmp{(i_comp-1)*2 + 1} = comps.gSac(i_comp).condition1;
            tmp{(i_comp-1)*2 + 2} = comps.gSac(i_comp).condition2;
        end
        gSacConditionFields = unique(tmp, 'stable');
    end

    % determine trials to include for gSac analysis
    % rule: include any trial that matches at least one condition AND is in g4f.analysisTrials (if present)
    if ~isempty(gSacConditionFields) && isfield(g4f, gSacConditionFields{1})
        nT = numel(g4f.(gSacConditionFields{1}));
        includedTrialsMat_gSac = false(nT, numel(gSacConditionFields));
        for i_col = 1:numel(gSacConditionFields)
            fld = gSacConditionFields{i_col};
            if isfield(g4f, fld)
                v = g4f.(fld);
                if numel(v) ~= nT
                    giveFeed(['warning: length mismatch for g4f field ' fld '. expected ' num2str(nT) ', got ' num2str(numel(v)) '. treating as empty.']);
                else
                    includedTrialsMat_gSac(:, i_col) = logical(v);
                end
            else
                giveFeed(['warning: condition field ' fld ' not found in g4f.']);
            end
        end
        gSacIncludedTrials = any(includedTrialsMat_gSac, 2);
    else
        % fallback: if we at least know how many trials from an event time vector, use all
        if isfield(g4f, 'targetOnTime')
            gSacIncludedTrials = true(numel(g4f.targetOnTime), 1);
            giveFeed('warning: gSacConditionFields empty; defaulting gSacIncludedTrials to all g4f trials.');
        else
            gSacIncludedTrials = false(0,1);
            giveFeed('warning: gSacConditionFields empty and no event time vector; gSacIncludedTrials is empty.');
        end
    end

    % intersect with analysisTrials if available
    if isfield(g4f, 'analysisTrials') && ~isempty(g4f.analysisTrials)
        if numel(g4f.analysisTrials) == numel(gSacIncludedTrials)
            gSacIncludedTrials = gSacIncludedTrials & logical(g4f.analysisTrials);
        else
            giveFeed('warning: g4f.analysisTrials length mismatch; not intersecting.');
        end
    end

    % create configuration for each gSac event
    for i_event = 1:numel(gSacEventLabels)
        cfg = struct();
        eventName = gSacEventLabels{i_event};
        cfg.eventName = eventName;

        if isfield(g4f, [eventName 'Time'])
            cfg.eventTimes = g4f.([eventName 'Time']);    % times for this event across all trials
        else
            giveFeed(['warning: ' eventName 'Time not found in g4f. setting to empty.']);
            cfg.eventTimes = [];
        end

        cfg.preTime         = gSacPreTimes(i_event);
        cfg.postTime        = gSacPostTimes(i_event);
        cfg.conditionFields = gSacConditionFields;        % e.g., 'isFaceTrial', 'isHighProbTrial', etc.
        cfg.includedTrials  = gSacIncludedTrials;         % logical vector for which trials to include
        cfg.analysisNeurons = data.analysisNeurons;       % selected neurons mask
        cfg.analysisTrials  = gSacIncludedTrials;         % alias for downstream fns

        % stash condition vectors for the included trials
        cfg.condVars = struct();
        if ~isempty(gSacIncludedTrials) && any(gSacIncludedTrials)
            for j_field = 1:numel(gSacConditionFields)
                fld = gSacConditionFields{j_field};
                if isfield(g4f, fld)
                    v = g4f.(fld);
                    if numel(v) == numel(gSacIncludedTrials)
                        cfg.condVars.(fld) = logical(v(:)) & gSacIncludedTrials(:);
                    else
                        giveFeed(['warning: length mismatch for g4f field ' fld '. skipping in condVars.']);
                    end
                end
            end
        end

        analysisConfig.(['sac_' eventName]) = cfg;        % store config, prefixed with 'sac_'
    end

    %% attention task (unchanged)
    attnEventLabels = {'fixOn','cueOn','stimOn','stimChg','reward'};
    attnPreTimes    = [-0.25, -0.10, -0.10, -0.10, -0.10];
    attnPostTimes   = [ 0.75,  0.90,  0.90,  0.90,  0.90];

    % collect unique condition fields from attn comparisons
    if isempty(comps.attn)
        attnConditionFields = {};
        giveFeed('warning: no attn comparisons defined. attn analysis will be limited.');
    else
        nAttnComps = numel(comps.attn);
        tmpA = cell(1, nAttnComps * 2);
        for i_comp = 1:nAttnComps
            tmpA{(i_comp-1)*2 + 1} = comps.attn(i_comp).condition1;
            tmpA{(i_comp-1)*2 + 2} = comps.attn(i_comp).condition2;
        end
        attnConditionFields = unique(tmpA, 'stable');
    end

    % determine trials to include for attn analysis
    if ~isempty(attnConditionFields) && isfield(attn, attnConditionFields{1})
        nT = numel(attn.(attnConditionFields{1}));
        includedTrialsMat_attn = false(nT, numel(attnConditionFields));
        for i_col = 1:numel(attnConditionFields)
            fld = attnConditionFields{i_col};
            if isfield(attn, fld)
                v = attn.(fld);
                if numel(v) ~= nT
                    giveFeed(['warning: length mismatch for attn field ' fld '. expected ' num2str(nT) ', got ' num2str(numel(v)) '. treating as empty.']);
                else
                    includedTrialsMat_attn(:, i_col) = logical(v);
                end
            else
                giveFeed(['warning: condition field ' fld ' not found in attn structure.']);
            end
        end
        attnIncludedTrials = any(includedTrialsMat_attn, 2);
    else
        if isfield(attn, 'stimOnTime')
            attnIncludedTrials = true(numel(attn.stimOnTime),1);
            giveFeed('warning: attnConditionFields empty; defaulting attnIncludedTrials to all attn trials.');
        else
            attnIncludedTrials = false(0,1);
            giveFeed('warning: attnConditionFields empty and no event time vector; attnIncludedTrials is empty.');
        end
    end

    % create configuration for each attn event
    for i_event = 1:numel(attnEventLabels)
        cfg = struct();
        eventName = attnEventLabels{i_event};
        cfg.eventName = eventName;

        if isfield(attn, [eventName 'Time'])
            cfg.eventTimes = attn.([eventName 'Time']);
        else
            giveFeed(['warning: ' eventName 'Time not found in attn. setting to empty.']);
            cfg.eventTimes = [];
        end

        cfg.preTime         = attnPreTimes(i_event);
        cfg.postTime        = attnPostTimes(i_event);
        cfg.conditionFields = attnConditionFields;
        cfg.includedTrials  = attnIncludedTrials;
        cfg.analysisNeurons = data.analysisNeurons;
        cfg.analysisTrials  = attnIncludedTrials;

        cfg.condVars = struct();
        if ~isempty(attnIncludedTrials) && any(attnIncludedTrials)
            for j_field = 1:numel(attnConditionFields)
                fld = attnConditionFields{j_field};
                if isfield(attn, fld)
                    v = attn.(fld);
                    if numel(v) == numel(attnIncludedTrials)
                        cfg.condVars.(fld) = logical(v(:)) & attnIncludedTrials(:);
                    else
                        giveFeed(['warning: length mismatch for attn field ' fld '. skipping in condVars.']);
                    end
                end
            end
        end

        analysisConfig.(['attn_' eventName]) = cfg;       % store config, prefixed with 'attn_'
    end

    % finalize
    data.analysisConfig = analysisConfig;
end