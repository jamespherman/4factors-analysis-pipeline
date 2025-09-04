function output_mat_file_path = local_computeArraysForAnalysisScOnly(matFileNamesSC_fullPaths_cell, giveFeed, oneDriveDir)
    % processes a single SC MAT file from Stage 1.5 to compute arrays for analysis,
    % selecting only "good" units based on cluster_KSLabel.tsv
    %
    % now supports both legacy 2-factor gSac and new 4-factor gSac (face/nonface/bullseye × sal × prob × reward-hemifield).
    % attention task is untouched.

    nFiles = numel(matFileNamesSC_fullPaths_cell);
    giveFeed('Stage 2 - Initializing parameters for array computation...');
    binWidth = 0.025; % seconds

    % prepare output directory
    neuronalDataAnalysisDir = fullfile(oneDriveDir, 'Neuronal Data Analysis');
    matDir2_output = fullfile(neuronalDataAnalysisDir, 'Newton Datasets');
    if ~exist(matDir2_output, 'dir')
        mkdir(matDir2_output);
        giveFeed(['Stage 2 - Created directory: ' matDir2_output]);
    end

    for i_file = 1:nFiles
        scFileFullPath = matFileNamesSC_fullPaths_cell{i_file};
        [currentDir, scFileNameOnly, scFileExtOnly] = fileparts(scFileFullPath);
        current_scFile_name = [scFileNameOnly scFileExtOnly];
        giveFeed(['Stage 2 - Processing file ' num2str(i_file) ' of ' num2str(nFiles) ': ' current_scFile_name]);

        % load stage 1.5 output
        giveFeed('Stage 2 - Loading data from Stage 1.5 output...');
        sc = load(scFileFullPath);

        % verify trialInfo exists
        if ~isfield(sc, 'trialInfo') || ~isfield(sc.trialInfo, 'nStim')
            error(['Stage 2 - ERROR: trialInfo or trialInfo.nStim missing in ' current_scFile_name]);
        end

        % clean up fields
        giveFeed('Stage 2 - Cleaning extraneous fields...');
        fieldsKeep = {'clusterCounts','clusterID','clusterTimes','eventTimes', ...
            'eventValuesTrials','fileDate','globalEventTimes','spikeClusters', ...
            'spikeTimes','trialInfo','nClusters','amplitudes','spikeTimesSamples', ...
            'initial_matFileName'};
        sc = rmfield(sc, setdiff(fieldnames(sc), fieldsKeep));

        % re-index clusters and compute counts/times
        giveFeed('Stage 2 - Re-indexing clusters and computing counts/times...');
        if isempty(sc.spikeClusters)
            sc.nClusters = 0;
            sc.clusterTimes = {};
        else
            [~,~,sc.spikeClusters] = unique(sc.spikeClusters);
            sc.nClusters = max(sc.spikeClusters);
            sc.clusterTimes = cell(1, sc.nClusters);
            for k = 1:sc.nClusters
                cid = sc.clusterCounts(k,1);
                sc.clusterTimes{k} = sc.spikeTimes(sc.spikeClusters == cid);
            end
        end

        % select only 'good' units from the same directory as the mat file
        giveFeed('Stage 2 - Selecting "good" units based on label file...');
        labelsFile = fullfile(currentDir, 'cluster_KSLabel.tsv');
        if exist(labelsFile, 'file')
            labelTable = readtable(labelsFile, 'FileType', 'text', 'Delimiter', '\t');
            goodUnits = labelTable{strcmp(labelTable{:,2}, 'good'), 1};
        else
            error(['Stage 2 - ERROR: label file missing: ' labelsFile]);
        end
        if sc.nClusters > 0
            qualityMask = false(sc.nClusters,1);
            validUnits = goodUnits(goodUnits >= 1 & goodUnits <= sc.nClusters);
            qualityMask(validUnits) = true;
        else
            qualityMask = false(0,1);
        end
        sc.analysisNeurons = qualityMask;
        giveFeed(['Stage 2 - Selected ' num2str(nnz(sc.analysisNeurons)) ' "good" units.']);

        % define event windows
        giveFeed('Stage 2 - Defining task event time windows...');
        sc.taskEventTimeWindows.sac_targetOn   = [0.05, 0.25];
        sc.taskEventTimeWindows.attn_stimOn    = [0.05, 0.25];
        sc.taskEventTimeWindows.sac_saccadeOn  = [-0.1, 0.05];

        % figure out which gSac flavor we have
        giveFeed('Stage 2 - Detecting gSac task flavor (2-factor vs 4-factor)...');
        is4F = false;
        try
            c = initCodes();
            if isfield(sc.trialInfo,'taskCode') && isfield(c,'uniqueTaskCode_gSac_4factors')
                is4F = any(sc.trialInfo.taskCode == c.uniqueTaskCode_gSac_4factors);
            end
        catch
            % if initCodes not on path, we’ll fall back later via targetOn presence; keep is4F=false
        end
        sc.is4F = is4F;

        % initialize task structures and side info
        giveFeed('Stage 2 - Initializing task structures...');
        % determine scSide using legacy gsac init (robust for both flavors)
        gSac_tmp = initializeGSacContrast(sc.trialInfo, sc.eventTimes);
        sc.scSide = local_defineScSide(sc, gSac_tmp, giveFeed);
        giveFeed(['Stage 2 - SC side: ' sc.scSide]);

        % re-init with scSide for downstream convenience
        if sc.is4F
            % 4-factor gSac (initialize4Factors internally calls decodeTaskCodes)
            g4f = initialize4Factors(sc.trialInfo, sc.eventTimes, sc.scSide);
            sc.g4f = g4f;
            gOut = g4f;  % unified handle passed to analysis config (helpers check sc.is4F)
        else
            % legacy 2-factor gSac
            gSac = initializeGSacContrast(sc.trialInfo, sc.eventTimes, sc.scSide);
            sc.gSac = gSac;
            gOut = gSac;
        end

        % attention task stays the same
        attn = initializeAttn(sc.trialInfo, sc.eventTimes, sc.eventValuesTrials, sc.scSide);

        % define comparisons and analysis config
        giveFeed('Stage 2 - Defining comparisons and analysis config...');
        % comparisons may differ by flavor; let the helper inspect sc.is4F and gOut
        sc.comparisons = local_defineComparisons();                 % helper should branch by flag
        sc = local_defineAnalysisConfig(gOut, attn, sc, giveFeed);        % helper should read sc.is4F and g4f/gSac

        % build PSTH data arrays
        giveFeed('Stage 2 - Building binned spike count arrays...');
        eventNames = fieldnames(sc.analysisConfig);
        for j = 1:numel(eventNames)
            sc = extractConditionMatrices(sc, eventNames{j}, binWidth);
        end

        % collect bin centers for plotting
        giveFeed('Stage 2 - Creating bin centers struct...');
        binCenters_struct = struct();
        for f = fieldnames(sc.analysisConfig)'
            ev = f{1};
            if isfield(sc.analysisConfig.(ev), 'binCenters') && ~isempty(sc.analysisConfig.(ev).binCenters)
                binCenters_struct.(ev) = sc.analysisConfig.(ev).binCenters;
            end
        end
        sc.binCenters = binCenters_struct;

        % save processed data
        output_mat_file_path = fullfile(matDir2_output, current_scFile_name);
        giveFeed(['Stage 2 - Saving processed data to: ' output_mat_file_path]);
        save(output_mat_file_path, '-struct', 'sc', '-v7.3');
        giveFeed(['Stage 2 - Done: ' current_scFile_name]);
    end
    giveFeed('--- Stage 2 - All files processed ---');
end
