clear
close all

% make in-line function to report timing / give feedback to user:
tic;
giveFeed = @(x)disp([num2str(toc) ' - ' x]);

% define mat file directory and save directory definition:
oneDriveDir = findOneDrive;
neuronalDataAnalysisDir = fullfile(oneDriveDir, 'Neuronal Data Analysis');
matDir = [neuronalDataAnalysisDir filesep '.mat analysis files' filesep];
saveFolder = [neuronalDataAnalysisDir filesep 'R01 2025 Figures\Individual Neuron Summary Plots'];

% list candidate .mat files:
fileCandidates = mydir(matDir);

% select the files we want - Feynman, recorded in August 2024:
matFileNamesSC = fileCandidates(...
    contains(fileCandidates, "Feynman") & ...
    contains(fileCandidates, "08") & ...
    contains(fileCandidates, "SC"));
matFileNamesSNc = fileCandidates(...
    contains(fileCandidates, "Feynman") & ...
    contains(fileCandidates, "08") & ...
    contains(fileCandidates, "SNc"));


% Determine the number of file pairs.
nFilePairs = numel(matFileNamesSC);

% give fedback
giveFeed('Analyzing sessions');

% Loop over file pairs.
for i = 1:nFilePairs
    % Feedback: processing file pair i.
    giveFeed(['Processing file pair ' num2str(i) ' of ' ...
        num2str(nFilePairs)]);

    % Define current SC and SNc filenames.
    scFile  = matFileNamesSC{i};
    sncFile = matFileNamesSNc{i};

    scData = load(fullfile(matDir, scFile));
    sncData = load(fullfile(matDir, sncFile));

    % extract arrays:
    plotIndividualNeuronPSTHs(scData, scFile, saveFolder);
    plotIndividualNeuronPSTHs(sncData, sncFile, saveFolder);

end

giveFeed('ALL DONE!');