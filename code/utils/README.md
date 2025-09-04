# Utility Functions

This directory contains robust, general-purpose utility functions that are intended to be called directly by the main analysis pipeline.

*   `alignAndBinSpikes.m`: Aligns spike times to given event times and bins them to create a histogram.
*   `arrayROC.m`: Computes the area under the Receiver Operating Characteristic (ROC) curve for binned spike count data.
*   `barStairsFill.m`: Plots two traces in a "stairs" style and fills the area between them with gray.
*   `calculate_baseline_fr.m`: Computes the average baseline firing rate for each neuron over a defined pre-trial period.
*   `calculate_waveform_metrics.m`: Analyzes a single mean waveform to extract features like the peak-to-trough duration.
*   `chi2proptest.m`: Performs a chi-square test to compare multiple proportions.
*   `code2string.m`: Converts numeric trial event codes into human-readable string labels.
*   `findOneDrive.m`: Finds the local path to the user's OneDrive directory, accommodating different operating systems.
*   `initCodes.m`: Initializes a "holy" set of numeric event codes used for synchronizing experiment events with electrophysiology recordings.
*   `mySubPlot.m`: Creates a subplot with precise, customizable control over grid position, margins, and spacing.
*   `mybarerr.m`: Creates a bar plot with customizable error bars and colors.
*   `outerLims.m`: Calculates the outermost x and y axis limits required to encompass a group of subplots.
*   `pdfSave.m`: Saves a specified figure to a PDF with WYSIWYG ("what you see is what you get") scaling.
*   `pngSave.m`: Saves a specified figure to a PNG with WYSIWYG ("what you see is what you get") scaling.
*   `richColors.m`: Provides a standardized palette of 16 aesthetically pleasing colors for plotting.
