# Data Dictionary for `session_data.mat`

This document provides a comprehensive description of the data structure contained within the `session_data.mat` file. This file is the final output of the data processing pipeline, created by `consolidate.consolidate_data.m`. It merges behavioral data (from `prep.prepare_behavioral_data.m`), spike data (from Kilosort), and mean waveform data (from `consolidate.extract_waveforms.m`).

The `session_data.mat` file contains a single top-level struct named `session_data`.

## `session_data` Structure

| Field | Type | Description |
|---|---|---|
| `trialInfo` | `table` | A table containing trial-by-trial information, merging behavioral parameters from PLDAPS with trial events from the neural data. Each row corresponds to one trial. |
| `eventTimes` | `struct` | A structure containing vectors of timestamps for key events in each trial. Each field corresponds to an event, and the vector length matches the number of trials. |
| `eventValuesTrials` | `cell` | A cell array where each cell contains a vector of the raw digital event codes (strobes) from the `.nev` file for the corresponding trial. |
| `spikes` | `struct` | A structure containing all spike-related data, including spike times, cluster assignments, and waveform information. |

---

## `session_data.trialInfo`

The `trialInfo` table is the primary source of behavioral data for each trial. The columns of this table are dynamically generated based on the specific PLDAPS task file used for the session. It includes a combination of data strobed to the neural recording system and parameters saved by the PLDAPS task code (`p.trVars` and `p.trData`).

Fields that are consistently scalar and numeric (or boolean) across all trials are stored as numeric arrays. Any trial where the value was empty (`[]`) in the raw PLDAPS data is represented as `NaN` in these arrays. Fields with inconsistent sizes (e.g., arrays of different lengths in different trials) or non-numeric data types are stored as `cell` arrays.

### Standard Fields (from `.nev` file)

These fields are derived from the `.nev` file and provide the basic structure of the trials.

| Field | Description |
|---|---|
| `trialNumber` | The sequential trial number from the `.nev` file. |
| `nevTrialStart` | The start time of the trial in neural recording clock samples. |
| `nevTrialEnd` | The end time of the trial in neural recording clock samples. |
| `nevTrialDuration` | The duration of the trial in seconds. |
| `pdsTrialStart` | The start time of the trial as recorded by PLDAPS, aligned to the neural recording clock. |
| `pdsTrialEnd` | The end time of the trial as recorded by PLDAPS, aligned to the neural recording clock. |
| `pdsTrialDuration` | The duration of the PLDAPS trial in seconds. |

### PLDAPS Fields (from `p.trVars` and `p.trData`)

**Note on Dynamic Fields**: The fields listed below are sourced from the PLDAPS `p` structure, which is saved on a per-trial basis. The `prepare_behavioral_data.m` script dynamically discovers all scalar and cell-based variables within the `p.trVars` and `p.trData` structures and adds them to the `trialInfo` table. Because different behavioral tasks save different variables, the list below should be considered a representative superset of possible fields, not an exhaustive list. Not all fields will be present in every `session_data.mat` file.

#### General & Control Variables

| Field | Description |
|---|---|
| `blockNumber` | The current block number. |
| `isVisSac` | Boolean flag, true if the trial was a visually-guided saccade trial. |
| `isCueChangeTrial` | Boolean flag, true if the trial was a cue change trial. |
| `isFoilChangeTrial` | Indicates if this is a foil change trial (1), no change (0), or if a foil is not present (-1). |
| `isNoChangeTrial` | Boolean flag, true if the trial was a no-change trial. |
| `isStimChangeTrial` | Boolean flag, true if it is a stimulus change trial. |
| `isStimChgNoDim` | Boolean, true if it is a stimulus change trial with no dimming. |
| `isContrastChangeTrial` | Boolean, true if it is a contrast change trial. |
| `passJoy` | If set to 1, simulates correct joystick behavior. Used for debugging. |
| `passEye` | If set to 1, simulates correct eye movement behavior. Used for debugging. |
| `repeat` | A flag that, if true, indicates the trial should be repeated. |
| `rwdJoyPR` | Determines reward condition (0 for joystick press, 1 for joystick release). |
| `trialEndState` | The numerical identifier of the state in which the trial ended. |
| `trialRepeatFlag` | A boolean flag that is true if the trial was an error and should be repeated. |
| `stimSeed` | Random seed for the stimulus parameters for the current trial. |
| `trialSeed` | Random seed for the trial parameters for the current trial. |

#### Stimulus, Geometry & Position Variables

| Field | Description |
|---|---|
| `fixDegX`, `fixDegY` | The X and Y coordinates of the fixation point in degrees of visual angle. |
| `fixLocRandX`, `fixLocRandY` | The range of random jitter for the fixation point location. |
| `fixWinWidthDeg`, `fixWinHeightDeg` | The width and height of the fixation window in degrees. |
| `targDegX`, `targDegY` | The X and Y coordinates of the target in degrees of visual angle. |
| `targTheta`, `targRadius` | The polar coordinates (angle and eccentricity) of the target. |
| `targWinWidthDeg`, `targWinHeightDeg` | The width and height of the target window in degrees. |
| `stimLoc1Elev`, `stimLoc1Ecc` | The elevation and eccentricity of the first stimulus location. |
| `stimLoc2Elev`, `stimLoc2Ecc`, etc. | The elevation and eccentricity for the other stimulus locations. |
| `stimElevs`, `stimEccs` | Vectors of elevations and eccentricities for all stimulus locations. |
| `stimOnList` | A numerical list of the stimuli that are on for the current trial. |
| `nStim` | The number of stimuli presented on the current trial. |
| `chgLoc` | The location of the stimulus change. |
| `cueLoc` | The location of the cue. |
| `changeDelta` | The magnitude of the stimulus change in the current trial. |
| `dimVals` | A vector of dimming values for each trial. |
| `lowDimVal`, `midDimVal`, `highDimVal` | The brightness levels for the fixation point after dimming. |
| `propVis` | The proportion of trials that are visually-guided. |
| `rwdSize` | The reward size for the current trial (e.g., 1=high, 2=low). |
| `stimType` | The type of stimulus (e.g., 1: Face, 2: Non-Face). |
| `salience` | The salience of the stimulus (e.g., 1 for high, 2 for low). |
| `reward` | The reward magnitude (e.g., 1 for high, 2 for low). |
| `preSacXY` | The average X,Y eye position in the 5ms before saccade onset. |
| `postSacXY` | The average X,Y eye position in the 5ms after saccade offset. |
| `peakVel` | The peak velocity of the primary saccade. |
| `SRT` | Saccadic Reaction Time, calculated as the time from the 'go' signal to saccade onset. |

#### Timing & Duration Variables

| Field | Description |
|---|---|
| `rewardDurationMs` | The duration of the reward delivery in milliseconds. |
| `rewardDelay` | The delay between a correct response and the delivery of the reward. |
| `timeoutdur` | The duration of the timeout period after an error. |
| `joyWaitDur` | The maximum time to wait for the joystick to be pressed at the start of a trial. |
| `fixWaitDur` | The maximum time to wait for the subject to acquire fixation. |
| `fixDur` | The required duration for which the subject must hold fixation. |
| `trialMax` | The maximum possible duration of a single trial. |
| `targHoldDuration` | The specific target hold duration for the current trial. |
| `targetFlashDuration` | For memory-guided trials, the duration the target is visible before disappearing. |
| `goLatencyMin`, `goLatencyMax` | The minimum and maximum allowed reaction time for initiating a saccade. |
| `timeTargOnset` | The specific time from fixation acquisition to target onset for the current trial. |
| `timeFixOffset` | The specific time from fixation acquisition to the 'go' signal for the current trial. |
| `stim2ChgIntvl` | The minimum time between stimulus onset and change. |
| `chgWinDur` | The time window during which a change can occur. |
| `joyMinLatency`, `joyMaxLatency` | The minimum and maximum acceptable joystick release latency. |

---

## `session_data.eventTimes`

The `eventTimes` structure contains vectors of timestamps for key trial events. Each field is a vector where the i-th element corresponds to the time of that event in the i-th trial. All timestamps are in seconds, aligned to the master neural recording clock.

### Standard Fields (from `.nev` file)

| Field | Description |
|---|---|
| `trialStart` | The timestamp of the trial's start. |
| `trialEnd` | The timestamp of the trial's end. |
| `fixOn` | Timestamp of when the fixation point appeared. |
| `fixAq` | Timestamp of when fixation was acquired. |
| `targOn` | Timestamp of when the saccade target appeared. |
| `targOff` | Timestamp of when the saccade target disappeared (in memory-guided trials). |
| `sacOn` | Timestamp of saccade onset. |
| `sacOff` | Timestamp of saccade offset. |
| `reward` | Timestamp of when the reward was delivered. |
| `outcomeOn` | For 'tokens' task trials, the corrected timestamp of when the token outcome was displayed. `NaN` for other tasks. |
| `brokeFix` | Timestamp of when a fixation break occurred. |
| `brokeJoy` | Timestamp of when a joystick break occurred. |

### PLDAPS Timing Fields (prefixed with `pds`)

**Note on Dynamic Fields**: These fields are dynamically discovered from the `p.trData.timing` structure in the PLDAPS data. The `prepare_behavioral_data.m` script identifies all scalar fields within this structure and adds them to `eventTimes`, prefixing the original field name with `pds`. For example, `p.trData.timing.targetReillum` becomes `eventTimes.pdsTargetReillum`. The list below is representative, but the exact fields present may vary between tasks.

| Field | Description |
|---|---|
| `pdsTrialStartPTB` | The trial start time as recorded by Psychtoolbox (`GetSecs`). |
| `pdsTrialStartDP` | The trial start time as recorded by the DataPixx. |
| `pdsTrialBegin` | Timestamp of when the trial logic began. |
| `pdsTrialEnd` | Timestamp of when the trial logic ended. |
| `pdsLastFrameTime`| The timestamp of the last screen flip. |
| `pdsFixOn` | Timestamp of when the fixation point appeared. |
| `pdsFixAq` | Timestamp of when fixation was acquired. |
| `pdsFixOff` | Timestamp of when the fixation point disappeared (the 'go' signal). |
| `pdsTargetOn` | Timestamp of when the saccade target appeared. |
| `pdsTargetOff` | Timestamp of when the saccade target disappeared. |
| `pdsTargetReillum`| Timestamp of when the target reappeared. |
| `pdsTargetAq` | Timestamp of when the target was acquired. |
| `pdsSaccadeOnset`| Timestamp of when the saccade was initiated. |
| `pdsSaccadeOffset`| Timestamp of when the saccade ended. |
| `pdsBrokeFix` | Timestamp of when a fixation break occurred. |
| `pdsReward` | Timestamp of when the reward was delivered. |
| `pdsTone` | Timestamp of when the auditory feedback was delivered. |
| `pdsJoyPress` | Time of joystick press. |
| `pdsJoyRelease` | Time of joystick release. |
| `pdsStimOn` | Timestamp of when the stimuli appeared. |
| `pdsStimOff` | Timestamp of when the stimuli disappeared. |
| `pdsCueOn` | Timestamp of when the cue appeared. |
| `pdsCueOff` | Timestamp of when the cue disappeared. |
| `pdsStimChg` | Timestamp of the stimulus change. |
| `pdsOutcomeOn` | Timestamp the token outcome was displayed. |

---

## Special Case: Timestamp Correction for `gSac_4factors` Task

In sessions that include the `gSac_4factors` task, a special correction procedure is applied during the `prep.prepare_behavioral_data` step to account for known timing inaccuracies in the data generated by this task.

### Correction Method

1.  **Identify Ground Truth**: The pipeline identifies all trials within the session that are *not* `gSac_4factors` trials to serve as a reliable "ground truth" for timestamp alignment.

2.  **Fit Linear Model**: It builds a linear regression model (`polyfit`) to map the PLDAPS-based timestamps from the reliable trials to the master Ripple-based timestamps.
    *   **Outlier Removal**: Before fitting the final model, an initial fit is performed to identify and remove outliers. Any data points where the residual (error) is more than 3 standard deviations from the mean are excluded from the final fit.

3.  **Predict and Correct Timestamps**: This linear model is then used to predict the correct, Ripple-based timestamps for all event times recorded within the `gSac_4factors` trials. The original, inaccurate timestamps for these trials are overwritten with the corrected values. This correction is applied to a specific set of event fields.

4.  **Null Unreliable Fields**: Certain event fields that are known to be particularly unreliable in this task are explicitly set to `NaN` for all `gSac_4factors` trials. These fields include `nonStart` and `blinkDuringSac`.

### Corrected Fields

The following fields in the `eventTimes` structure are re-calculated based on the linear model for all `gSac_4factors` trials:

*   `CUE_ON`
*   `fixAq`
*   `fixBreak`
*   `fixOff`
*   `fixOn`
*   `joyPress`
*   `lowTone`
*   `reward`
*   `REWARD_GIVEN`
*   `saccadeOffset`
*   `saccadeOnset`
*   `targetAq`
*   `targetOff`
*   `targetOn`
*   `targetReillum`
*   `trialEnd`
*   `TRIAL_END`
*   `trialBegin`

### Diagnostic Plot

As part of this process, a diagnostic plot is generated to allow for visual inspection of the quality of the linear fit. This plot, named `{unique_id}_timestamp_fit.pdf`, is saved in the `diagnostics` subfolder within the session's processed data directory (e.g., `.../{unique_id}/diagnostics/`). It shows a scatter plot of the ground-truth timestamps (after outlier removal) and the best-fit line, along with the R-squared value of the fit.

---

## Special Case: Timestamp Correction for 'tokens' Task Outcome

For trials belonging to the 'tokens' task, the general session-wide timestamp correction model (as described in the `gSac_4factors` section) is used to calculate a precise, corrected timestamp for when the trial's outcome (i.e., the token) was displayed.

### Correction Method

1.  **Identify 'tokens' Trials**: The pipeline identifies all trials where the `taskCode` corresponds to the 'tokens' task.
2.  **Apply Linear Model**: For each of these trials, the pipeline takes the PLDAPS-based timestamp for the outcome display (`pdsOutcomeOn`) and uses the session's linear correction model to translate it into the master Ripple clock's time base.
3.  **Store Corrected Timestamp**: The resulting corrected timestamp is stored in a new, dedicated field in the `eventTimes` structure called `outcomeOn`. If a trial is not a 'tokens' task trial, its value in `eventTimes.outcomeOn` will be `NaN`.

---

## `session_data.spikes`

The `spikes` structure contains all data related to neural spikes, sourced from the output of Kilosort and subsequent analysis.

| Field | Type | Description |
|---|---|---|
| `times` | `Nx1 double` | An array of spike times in seconds, aligned to the master neural recording clock. `N` is the total number of spikes in the session. |
| `clusters` | `Nx1 uint32` | An array of cluster IDs corresponding to each spike in `spikes.times`. |
| `cluster_info` | `table` | A table containing metrics and information about each sorted cluster, imported from Kilosort's `cluster_info.tsv` file. |
| `wfMeans` | `cell` | A cell array where each cell contains the mean waveform for a cluster. The waveform is a `C x T` matrix, where `C` is the number of channels and `T` is the number of time samples around the spike. |
| `wfStds` | `cell` | A cell array where each cell contains the standard deviation of the waveforms for a cluster. The dimensions are the same as for `wfMeans`. |

### `spikes.cluster_info` Table

This table provides detailed information about each cluster identified by Kilosort. The exact columns can vary, but they typically include:

| Column | Description |
|---|---|
| `cluster_id` | The unique ID for the cluster. |
| `amplitude` | The average amplitude of the spikes in this cluster. |
| `contam_pct` | An estimate of the contamination percentage for this cluster. |
| `group` | The quality label assigned to the cluster (e.g., 'good', 'mua' for multi-unit activity, 'noise'). |
| `ch` | The primary channel on which the cluster's spikes were detected. |
| `depth` | The estimated depth of the cluster along the probe. |
| `fr` | The firing rate of the cluster in spikes per second. |
| `n_spikes` | The total number of spikes in the cluster. |
| `sh` | A measure of the refractory period violation rate. |
| `KSLabel` | The label assigned by the Kilosort algorithm itself (e.g., 'good' or 'mua'). |
