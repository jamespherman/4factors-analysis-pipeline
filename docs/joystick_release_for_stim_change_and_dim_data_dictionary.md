# Data Dictionary for 'joystick_release_for_stim_change_and_dim'

This document describes the data structures saved during the 'joystick_release_for_stim_change_and_dim' task, with a focus on the variant defined by the 'joystick_release_for_orient_change_and_dim_learn_cue_settings.m' settings file.

The data is saved in a main structure `p`. The entire `p` structure is saved once at the beginning of the session in a file named `p.mat`. For each trial, a separate file named `trialXXXX.mat` is saved, containing the `trVars`, `trData`, `status`, and `init` structures for that trial.

## The `p` Structure (Session-Level Data)

The `p` structure contains all the parameters and data for the experimental session. The fields described below are those that are saved once at the beginning of the session.

### `p.init`

This structure contains initialization parameters that are set once at the beginning of the session.

| Field | Description |
|---|---|
| `pcName` | The hostname of the computer running the experiment. |
| `rigConfigFile` | The path to the rig configuration file, which contains subject- and rig-specific details. |
| `exptType` | A string identifying the specific experiment being run. For this task variant, it is 'joystick_release_for_stim_dim_and_orient_change_learn_cue_multi'. This string determines which trial structure table is used. |
| `taskName` | The name of the task, which is 'joystick_release_for_stim_change_and_dim'. |
| `taskType` | A numerical index for the task type (value is 1). |
| `pldapsFolder` | The path to the main PLDAPS directory. |
| `protocol_title` | The banner text to identify the experimental protocol. |
| `date` | The date of the session in 'yyyymmdd' format. |
| `time` | The time of the session in 'HHMM' format. |
| `date_1yyyy`, `date_1mmdd`, `time_1hhmm` | Numerical representations of the date and time, with a '1' prepended to avoid issues with leading zeros. |
| `useDataPixxBool` | A boolean indicating whether the DataPixx/ViewPixx is being used. |
| `outputFolder` | The folder where the output files are saved. |
| `figureFolder` | The folder where figures are saved. |
| `sessionId` | A unique identifier for the session, combining the date, time, and task name. |
| `sessionFolder` | The folder where the data for the current session is saved. |
| `taskFiles` | A structure containing the filenames for the init, next, run, and finish scripts for the task. |
| `taskActions` | A cell array of strings containing the names of action M-files to be used in the task. |
| `trDataInitList` | A cell array that defines the variables within `p.trData` and their initial values at the start of each trial. |
| `nTrDataListRows` | The number of rows in `p.init.trDataInitList`. |
| `strobeList` | A cell array of variable names to be strobed at the end of each trial. |
| `trialsArray` | A matrix that defines the parameters for each trial in a block. Each row represents a trial, and each column represents a parameter. The columns are defined in `p.init.trialArrayColumnNames`. |
| `trialArrayColumnNames` | A cell array of strings that are the names of the columns in `p.init.trialsArray`. |
| `trialsTable` | The table of trial conditions used to generate `p.init.trialsArray`. |
| `blockLength` | The number of trials in a block. |
| `codes` | A structure containing the numerical codes for various trial events that are strobed to the ephys recording system. |
| `strb` | An object of the `pds.classyStrobe` class, used for strobing values to the ephys system. |

### `p.audio`

This structure contains parameters related to audio feedback.

| Field | Description |
|---|---|
| `audsplfq` | The audio playback sampling rate for the DataPixx (48000 Hz). |
| `Hitfq` | The frequency of the "hit" tone (600 Hz). |
| `Missfq` | The frequency of the "low" (miss) tone (100 Hz). |
| `auddur` | The duration of the tones in samples (4800 samples). |
| `lineOutLevel` | The audio level for the DataPixx line out (0 to 1). |
| `pcPlayback` | A boolean indicating whether to use the PC for audio playback. |

### `p.draw`

This structure contains parameters related to drawing visual elements on the screen.

| Field | Description |
|---|---|
| `ringThickDeg` | The thickness of the cue ring in degrees. |
| `ringRadDeg` | The radius of the cue ring in degrees. |
| `eyePosWidth` | The width of the eye position indicator in pixels. |
| `fixPointWidth` | The line width of the fixation point in pixels. |
| `fixPointRadius` | The radius of the fixation point in pixels. |
| `fixWinPenPre` | The pen width for the fixation window before a change event. |
| `fixWinPenPost` | The pen width for the fixation window after a change event. |
| `fixWinPenDraw` | The current pen width for the fixation window, which is set to either `fixWinPenPre` or `fixWinPenPost` during the trial. |
| `gridSpacing` | The spacing of the grid on the experimenter's display in degrees. |
| `gridW` | The grid spacing in degrees. |
| `joyRect` | The rectangle defining the position of the joystick indicator on the experimenter's display. |
| `cursorW` | The width of the cursor in pixels. |
| `clutIdx` | A structure containing integer indices for the Color Look-Up Table (CLUT). Each field name is a descriptive string for a color, and the value is the row index in the CLUT. |
| `color` | A structure that holds the current CLUT index for various visual elements (e.g., `p.draw.color.background`, `p.draw.color.fix`). The values of these fields are updated during the trial to change the colors of the elements. |
| `cueRingRect` | The rectangle defining the position and size of the cue ring in pixels. |
| `ringThickPix` | The thickness of the cue ring in pixels. |
| `myCLUTs` | The Color Look-Up Tables for each frame of the stimulus animation. |
| `stimTex` | A cell array of textures for the stimuli. |
| `middleXY` | The pixel coordinates of the center of the screen. |
| `gridXY` | The coordinates for the grid lines on the experimenter display. |
| `window` | The handle for the Psychtoolbox window. |
| `cueArcAngles` | The start and end angles for the colored and grey portions of the cue arc. |
| `cueArcProp` | The proportion of the cue arc that is colored, based on the probability of a change at the cued location. |
| `fixPointPix` | The pixel coordinates of the fixation point. |
| `fixPointRect` | The rectangle defining the fixation point. |

### `p.state`

This structure defines the numerical codes for the different states of the trial state machine.

| State | Code |
|---|---|
| `trialBegun` | 1 |
| `waitForJoy` | 2 |
| `showFix` | 3 |
| `dontMove` | 4 |
| `makeDecision` | 5 |
| `fixBreak` | 11 |
| `joyBreak` | 12 |
| `nonStart` | 13 |
| `hit` | 21 |
| `cr` | 22 |
| `miss` | 23 |
| `foilFa` | 24 |
| `fa` | 25 |

### `p.status`

This structure contains variables that track the status of the experiment across trials. It is updated at the end of each trial.

| Field | Description |
|---|---|
| `iTrial` | The current trial number. |
| `iGoodTrial` | The number of "good" trials (hits, misses, correct rejects, foil FAs). |
| `trialsLeftInBlock` | The number of trials remaining in the current block. |
| `blockNumber` | The current block number. |
| `fixDurReq` | The required fixation duration for the last trial. |
| `hr1stim`, `hr2stim`, `hr3stim`, `hr4stim` | Hit rates for trials with 1, 2, 3, and 4 stimuli, respectively. |
| `hc1stim`, `hc2stim`, `hc3stim`, `hc4stim` | Hit counts for trials with 1, 2, 3, and 4 stimuli, respectively. |
| `tc1stim`, `tc2stim`, `tc3stim`, `tc4stim` | Total trial counts for trials with 1, 2, 3, and 4 stimuli, respectively. |
| `hr1Loc1`, `cr1Loc1`, `hr1Loc2`, `cr1Loc2`, etc. | Hit rates and correct reject rates for different stimulus configurations. |
| `hc1Loc1`, `crc1Loc1`, `hc1Loc2`, `crc1Loc2`, etc. | Hit counts and correct reject counts for different stimulus configurations. |
| `cue1CtLoc1`, `foil1CtLoc1`, `cue1CtLoc2`, `foil1CtLoc2`, etc. | Counts of cue and foil change trials for different stimulus configurations. |
| `totalHits` | Total number of hits. |
| `totalMisses` | Total number of misses. |
| `totalChangeFalseAlarms` | Total number of false alarms on change trials. |
| `totalNoChangeFalseAlarms` | Total number of false alarms on no-change trials. |
| `totalCorrectRejects` | Total number of correct rejects. |
| `missedFrames` | The number of missed frames reported by Psychtoolbox. |
| `freeRwdRand` | The random number drawn to determine if a free reward should be given. |
| `freeRwdTotal` | The total number of free rewards delivered. |
| `freeRwdLast` | The trial number of the last free reward. |
| `trialsArrayRowsPossible` | A logical vector indicating which rows of the `trialsArray` are still available to be run in the current block. |
| `freeRewardsAvailable` | A logical vector indicating which trials in the block are designated to have a free reward. |
| `trialEndStates` | A vector containing the end state of each trial. |
| `reactionTimes` | A vector of reaction times for each trial. |
| `dimVals` | A vector of dimming values for each trial. |
| `changeDelta` | The magnitude of the stimulus change in the current trial. |
| `chgLoc` | The location of the stimulus change. |
| `cueLoc` | The location of the cue. |
| `nStim` | The number of stimuli on the current trial. |

### `p.trVars`

This structure contains the parameters for the current trial. Its values are set at the beginning of each trial in `_next.m`, primarily by copying from `p.trVarsGuiComm` and then being modified by `nextParams.m`.

| Field | Description |
|---|---|
| `passJoy` | If set to 1, simulates a correct joystick response. Used for debugging. |
| `passEye` | If set to 1, simulates correct eye fixation. Used for debugging. |
| `blockNumber` | The current block number. |
| `repeat` | If true, the trial will be repeated. |
| `rwdJoyPR` | Determines reward condition (0 for joystick press, 1 for joystick release). |
| `isCueChangeTrial` | Boolean indicating if the trial is a cue change trial. |
| `isFoilChangeTrial` | Indicates if this is a foil change trial (1), no change (0), or if a foil is not present (-1). |
| `isNoChangeTrial` | Boolean indicating if the trial is a no-change trial. Set in `nextParams.m`. |
| `finish` | The maximum number of trials to run. |
| `filesufix` | A suffix for the saved file. |
| `joyVolt` | The current joystick voltage. |
| `eyeDegX`, `eyeDegY` | The current eye position in degrees. |
| `eyePixX`, `eyePixY` | The current eye position in pixels. |
| `propHueChgOnly` | The proportion of change trials where only the hue of the peripheral stimulus changes (no dimming). |
| `isStimChangeTrial` | Boolean, true if it is a stimulus change trial. Set in `nextParams.m`. |
| `chgAndDimOnMultiOnly` | Boolean, if true, change+dim trials only occur on multi-stimulus trials. |
| `stimLoc1Elev`, `stimLoc1Ecc` | The elevation and eccentricity of the first stimulus location. |
| `stimLoc2Elev`, `stimLoc2Ecc`, etc. | The elevation and eccentricity for the other stimulus locations. Can be used to override the default circular arrangement. |
| `fixDegX`, `fixDegY` | The X and Y coordinates of the fixation point in degrees. |
| `fixLocRandX`, `fixLocRandY` | The range of random jitter for the fixation point location. |
| `lowDimVal`, `midDimVal`, `highDimVal` | The brightness levels for the fixation point after dimming, relative to the background. |
| `speedInit`, `ctrstInit`, `orientInit`, `freqInit`, `satInit`, `lumInit`, `hueInit` | Initial values for the various stimulus features. |
| `orientVar`, `hueVar`, `lumVar`, `satVar` | The variance for stimulus features that can be varied. |
| `speedDelta`, `contDelta`, `orientDelta`, `freqDelta`, `satDelta`, `lumDelta`, `hueDelta` | The magnitude of the change for each stimulus feature on change trials. |
| `stimRadius` | The radius of the stimulus aperture in degrees. |
| `boxSizePix` | The size of the "checks" in the checkerboard stimulus in pixels. |
| `boxLifetime` | The lifetime of the "checks" in frames. |
| `nPatches` | The number of stimulus patches. |
| `nEpochs` | The number of stimulus epochs (pre-change and post-change). |
| `rewardDurationMs` | The duration of the reward in milliseconds. |
| `rewardDurationMsSmall` | The duration of a small reward in milliseconds. |
| `fix2CueIntvlMin`, `fix2CueIntvlWin` | The minimum and window duration for the interval between fixation and cue onset. |
| `fix2CueIntvl` | The actual interval between fixation and cue onset for the current trial, randomly drawn from the min/win range. |
| `cueDur` | The duration of the cue presentation. |
| `cue2StimIntvlMin`, `cue2StimIntvlWin` | The minimum and window duration for the interval between cue offset and stimulus onset. |
| `cue2StimIntvl` | The actual interval between cue offset and stimulus onset for the current trial. |
| `stim2ChgIntvl` | The minimum time between stimulus onset and change. |
| `chgWinDur` | The time window during which a change can occur. |
| `rewardDelay` | The delay between a correct response (hit) and reward delivery. |
| `joyMinLatency`, `joyMaxLatency` | The minimum and maximum acceptable joystick release latency. |
| `timeoutAfterFa`, `timeoutAfterFoilFa`, `timeoutAfterMiss`, `timeoutAfterFixBreak` | The timeout durations for different trial outcomes. |
| `joyWaitDur` | The maximum time to wait for a joystick press at the start of a trial. |
| `fixWaitDur` | The maximum time to wait for fixation acquisition. |
| `freeDur` | The time before the start of the joystick press check. |
| `trialMax` | The maximum duration of a trial. |
| `joyReleaseWaitDur` | The time to wait after trial end to start the end-of-trial flicker if the joystick is not released. |
| `stimFrameIdx` | The current frame index for the stimulus animation. |
| `flipIdx` | The index of the current screen flip. |
| `postRewardDurMin`, `postRewardDurMax` | The minimum and maximum duration to wait after reward delivery before ending the trial. |
| `useQuest` | Boolean, if true, QUEST is used to determine stimulus parameters. |
| `numTrialsForPerfCalc` | The number of recent trials to use for performance calculation. |
| `freeRewardProbability` | The probability of a free reward between trials. |
| `freeRewardFlag` | A boolean to enable or disable free rewards. |
| `connectRipple` | Boolean, if true, connect to the Ripple system. |
| `rippleChanSelect` | The selected Ripple channel. |
| `useOnlineSort` | Boolean, if true, use online sorted spike times from Trellis. |
| `psthBinWidth`, `fixOnPsthMinTime`, etc. | Parameters for PSTH plotting. |
| `currentState` | The current state of the trial state machine. |
| `exitWhileLoop` | A boolean that controls the main trial loop. |
| `cueIsOn` | A boolean indicating if the cue is currently displayed. |
| `stimIsOn` | A boolean indicating if the stimuli are currently displayed. |
| `fixWinWidthDeg`, `fixWinHeightDeg` | The width and height of the fixation window in degrees. |
| `fixPointRadPix`, `fixPointLinePix` | The radius and line width of the fixation point in pixels. |
| `useCellsForDraw` | A boolean for drawing options. |
| `wantEndFlicker` | A boolean to enable the end-of-trial screen flicker. |
| `wantOnlinePlots` | A boolean to enable online plotting. |
| `fixColorIndex` | The color index for the fixation point. |
| `postFlip` | A structure used to log the timing of events that occur immediately after a screen flip. |
| `optoStimDurSec`, `optoPulseDurSec`, `otoPulseAmpVolts`, `optoIpiSec`, `isOptoStimTrial` | Parameters for optogenetic stimulation. |
| `currentTrialsArrayRow` | The row index of `p.init.trialsArray` for the current trial. |
| `stimSeed`, `trialSeed` | Random seeds for the stimulus and trial parameters. |
| `isStimChgNoDim` | Boolean, true if it is a stimulus change trial with no dimming. |
| `stim1On`, `stim2On`, etc. | Booleans indicating which stimuli are on for the current trial. |
| `stimOnList` | A numerical list of the stimuli that are on. |
| `isContrastChangeTrial` | Boolean, true if it is a contrast change trial. |
| `stimElevs`, `stimEccs` | Vectors of elevations and eccentricities for all stimulus locations. |
| `stimLocCart`, `stimLocCartPix` | The Cartesian coordinates of the stimulus locations in degrees and pixels. |
| `stimRects` | The rectangles defining the stimulus patches. |
| `fix2StimOnIntvl` | The interval between fixation acquisition and stimulus onset. |
| `stimChangeTime` | The time of the stimulus change (or pseudo-change) relative to fixation acquisition. |
| `joyMaxLatencyAfterChange` | The latest acceptable joystick release time after a change. |
| `hitRwdTime`, `corrRejRwdTime` | The time of reward delivery for hits and correct rejects. |
| `fix2StimOffIntvl` | The interval between fixation acquisition and stimulus offset. |
| `stimDur` | The maximum possible stimulus duration. |
| `stimFrames` | The total number of frames for the stimulus animation. |
| `rewardScheduleDur` | The duration of the reward schedule. |
| `postRewardDuration` | The duration to wait after reward delivery before ending the trial. |

### `p.trData`

This structure contains the data collected during the current trial.

| Field | Description |
|---|---|
| `eyeX`, `eyeY`, `eyeP`, `eyeT` | These fields are intended to store eye data, but in the current implementation, they are not populated. Eye data is stored in `p.trData.onlineEyeX` and `p.trData.onlineEyeY`. |
| `joyV` | Stores the joystick voltage values. Not currently populated in the provided code. |
| `dInValues`, `dInTimes` | Store the values and times of digital input events. Not currently populated in the provided code. |
| `spikeTimes` | Stores the timestamps of spikes from the ephys system. |
| `eventTimes`, `eventValues` | Store the times and values of strobed events. |
| `onlineEyeX`, `onlineEyeY` | Store the X and Y eye position in degrees, sampled on every frame of the trial. |
| `spikeClusters` | Stores the cluster ID for each spike. |
| `timing` | A sub-structure containing the timestamps of various trial events, relative to the start of the trial. See the **Timing Information** section below for a detailed explanation. |
| `trialEndState` | The numerical code of the state in which the trial ended. |
| `dimVal` | Set to 1 for change trials and 0 for no-change trials. This is a legacy variable and may not be used in analysis. |
| `reactionTime` | The reaction time, calculated as `joyRelease` - `stimChangeTime`. |
| `fixHoldReqMet` | The time the required fixation duration was met. |

### Timing Information (`p.trData.timing`)

This substructure contains timestamps for all critical events within the trial. Before each trial, all timing fields are initialized to -1. During the trial, they are populated with timestamps.

**All timestamps are measured in seconds relative to `p.trData.timing.trialStartPTB`**. This reference time is captured at the very beginning of the trial's `_run.m` function using `pds.getTimes`, which calls Psychtoolbox's `GetSecs`. Subsequent event times are calculated as `timeNow = GetSecs - p.trData.timing.trialStartPTB;`.

| Field | Description |
|---|---|
| `trialStartPTB` | The trial start time as recorded by Psychtoolbox (`GetSecs`). This is the master reference time for all other timestamps in this structure. |
| `trialStartDP` | The trial start time as recorded by the DataPixx. Can be used for synchronization. |
| `trialBegin` | The time the trial began (state 1). |
| `joyPress` | The time the joystick was pressed. |
| `fixOn` | The time the fixation point appeared. |
| `fixAq` | The time fixation was acquired. |
| `stimOn` | The time the stimuli appeared. |
| `stimOff` | The time the stimuli disappeared. |
| `cueOn` | The time the cue appeared. |
| `cueOff` | The time the cue disappeared. |
| `stimChg` | The time of the stimulus change. |
| `noChg` | The time of the pseudo-change on a no-change trial. |
| `brokeFix` | The time of a fixation break. |
| `brokeJoy` | The time of a joystick break. |
| `reward` | The time of reward delivery. |
| `tone` | The time of audio feedback. |
| `joyRelease` | The time the joystick was released. |
| `freeReward` | The time a free reward was delivered. |
| `flipTime` | A vector of timestamps for each screen flip. |
| `lastFrameTime` | The time of the last screen flip. |
| `optoStim` | The time of optogenetic stimulation onset. |
| `optoStimSham` | The time of sham optogenetic stimulation. |

### `p.stim`

This structure contains parameters and data related to the stimuli.

| Field | Description |
|---|---|
| `featureValueNames` | A cell array of strings with the names of the stimulus features that can be varied. |
| `nFeatures` | The number of features in `p.stim.featureValueNames`. |
| `stimChgIdx` | The index of the stimulus that changes on the current trial. |
| `cueLoc` | The index of the cued stimulus location. |
| `nStim` | The number of stimuli presented on the current trial. |
| `primStim` | The index of the "primary" stimulus, used to ensure that stimuli have different starting values for their features. |
| `speedArray`, `ctrstArray`, `orientArray`, `freqArray`, `satArray`, `lumArray`, `hueArray` | These are arrays that define the values of the stimulus features for each stimulus patch (rows) and each epoch (columns). The values are determined by the initial values in `p.trVars`, the deltas for the current trial, and the `trialsArray`. |
| `orientVarArray`, `hueVarArray`, `lumVarArray`, `satVarArray` | Arrays defining the variance for the corresponding stimulus features. |
| `patchDiamPix` | The diameter of the stimulus patch in pixels. |
| `patchDiamBox` | The diameter of the stimulus patch in "boxes" or "checks". |
| `epochFrames` | A vector containing the number of frames in each epoch. |
| `chgFrames` | A cumulative sum of the epoch durations in frames, used to determine when changes occur. |
| `funs` | A structure containing anonymous functions for various calculations. |
| `X`, `Y` | 3D arrays containing the X and Y indices for each box in each frame of the stimulus animation. |
| `nBoxTot` | The total number of boxes across all patches. |
| `uci` | A 3D array containing a unique numerical index for each unique color value needed for the stimulus animation. |
| `upi` | A 3D array containing the unique phase index for the gabor phase values. |
| `ucir` | A reshaped version of `p.stim.uci`. |
| `tempR`, `tempG`, `tempB` | Temporary vectors used to store the RGB color values for the stimuli. |
| `colorRowSelector` | A logical vector used to select the correct color values from the temporary color vectors. |
| `stimArray` | A cell array where each cell contains the image matrix for a stimulus patch. |
