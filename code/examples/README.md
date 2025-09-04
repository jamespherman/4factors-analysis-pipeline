# Example Scripts & Functions

This directory contains example scripts and functions that serve as conceptual blueprints, logical references, and stylistic patterns.

**IMPORTANT**: The code in this directory should not be called directly by the main analysis pipeline. Its purpose is to provide a reference for developing new, robust analysis and plotting functions.

*   `local_computeArraysForAnalysisScOnly.m`: Processes a single SC MAT file from Stage 1.5 to compute arrays for analysis, selecting only "good" units based on `cluster_KSLabel.tsv`.
*   `local_computeBaselineVsComparison.m`: Computes an ROC analysis to compare neural activity between a baseline period and a post-event period for contralateral and ipsilateral conditions.
*   `local_computeDiscriminability.m`: Computes both single-neuron (ROC) and population-level (SVM) discriminability between two experimental conditions over time.
*   `local_computeFiringRate.m`: Computes the average firing rate for a set of neurons over a defined time window.
*   `local_computePopulationDiscriminability.m`: Computes population-level discriminability between two conditions for each time bin using a support vector machine (SVM) classifier.
*   `local_computeWindowedROC.m`: Computes an ROC analysis for pre-defined comparisons on neural data that has been averaged within a specific time window for a given task event.
*   `local_defineAnalysisConfig.m`: Defines the configuration for various analysis events, specifying time windows, conditions, and trials to be included.
*   `local_defineComparisons.m`: Defines the sets of experimental conditions to be compared for both the gSac and attention tasks.
*   `local_defineScSide.m`: Determines the recording side of the superior colliculus (left or right) by comparing neural responses to contralateral and ipsilateral stimuli.
*   `local_selectNeurons.m`: Selects task-modulated neurons by comparing their firing rates across different, behaviorally-relevant time epochs (e.g., baseline, visual, saccade).
*   `makePSTHGrid.m`: Generates a grid of PSTH subplots for different groups of neurons, with options for custom colors and multiple event types.
*   `performCCA.m`: Performs a Canonical Correlation Analysis (CCA) to identify shared patterns of activity between two simultaneously recorded neural populations (SC and SNc).
*   `plotIndividualNeuronPSTHs.m`: Generates and saves a tiled figure for each individual neuron, showing its PSTHs for various task events and experimental comparisons.
*   `plotIndividualNeuronPSTHs_script.m`: Loads paired SC and SNc data files and calls `plotIndividualNeuronPSTHs` for each.
*   `plotPopulationPreferenceScOnly.m`: Creates a complex, multi-panel figure showing the proportion of SC neurons that prefer contralateral vs. ipsilateral stimuli across different task conditions and events.
*   `unitQualityAssesment_script.m`: Assesses the quality of single-unit recordings based on refractory period violations and mean firing rate.
