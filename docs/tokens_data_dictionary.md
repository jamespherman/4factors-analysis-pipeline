# Data Dictionary for 'tokens' Task

This document provides a comprehensive explanation of the fields and subfields of the data structures saved on each trial of the 'tokens' task. The primary data structure is `p`, which is saved once per session. This document will detail all sub-structures within `p`.

## `p.init`

This structure contains initialization parameters that are set only once at the very beginning of an experimental session.

| Field | Description |
|---|---|
| `pcName` | The hostname of the computer running the experiment. Used to select the correct rig configuration file. |
| `rigConfigFile` | The path to the rig-specific configuration file (e.g., `rigConfig_rig1.m`). This file contains hardware details like screen distance. |
| `useDataPixxBool` | A boolean flag that should always be `true`, indicating that the DataPixx I/O box is in use. |
| `taskName` | A string identifying the name of the task, which is `'tokens'`. |
| `taskType` | A numerical index for the task type (value is `1`). The documentation notes this is "poorly defined". |
| `pldapsFolder` | The root directory of the PLDAPS-based task. |
| `protocol_title` | A string used for the banner text in the GUI to identify the protocol (e.g., `'tokens_task'`). |
| `date` | The date the session was started, in `'yyyymmdd'` format. |
| `time` | The time the session was started, in `'HHMM'` format. |
| `date_1yyyy` | The year of the session, prepended with a '1' to avoid issues with leading zeros when converting to a double. |
| `date_1mmdd` | The month and day of the session, prepended with a '1'. |
| `time_1hhmm` | The hour and minute of the session, prepended with a '1'. |
| `outputFolder` | The path to the main output directory where session data is saved. |
| `figureFolder` | The path to the directory where figures generated during the session are saved. |
| `sessionId` | A unique identifier for the session, combining the date, time, and task name (e.g., `'20230525_t1430_tokens'`). |
| `sessionFolder` | The full path to the directory for the current session's data. |
| `taskFiles` | A structure containing the filenames for the core task functions (`init`, `next`, `run`, `finish`). |
| `taskActions` | A cell array of strings listing user-defined "action" M-files that can be called from the GUI. |
| `trialsPerCondition` | The number of times each unique trial condition is repeated within a block. |
| `exptType` | A string indicating the specific version of the experiment being run (e.g., `'tokens_AV'`). |
| `trDataInitList` | A cell array that defines the initial values for variables within `p.trData`. This is used to reset trial-specific data at the beginning of each new trial. |
| `nTrDataListRows` | The number of rows in `p.init.trDataInitList`, stored for efficient looping. |
| `strobeList` | A cell array defining a list of variables whose values are to be strobed (sent as event markers to the ephys system) at the end of each trial. Each row contains the variable to be strobed and a human-readable name for it. |

## `p.rig`

This structure contains parameters related to the specific hardware configuration of the experimental rig.

| Field | Description |
|---|---|
| `guiStatVals` | A cell array of strings listing status variables (`p.status`) to be displayed in the GUI. |
| `guiVars` | A cell array of strings listing trial variables (`p.trVarsInit`) that can be modified from the GUI. |
| `dp` | A sub-structure containing settings for the DataPixx. |
| `dp.useDataPixxBool` | A boolean indicating if the DataPixx is being used. |
| `dp.adcRate` | The sampling rate for the Analog-to-Digital Converter (ADC) in Hz (e.g., 1000). |
| `dp.maxDurADC` | The maximum duration in seconds to pre-allocate for the ADC buffer (e.g., 15). |
| `dp.adcBuffAddr` | The memory buffer address for the ADC on the DataPixx. |
| `dp.dacRate` | The sampling rate for the Digital-to-Analog Converter (DAC) in Hz (e.g., 1000). |
| `dp.dacPadDur` | A padding duration in seconds for the DAC signal. |
| `dp.dacBuffAddr` | The memory buffer base address for the DAC on the DataPixx. |
| `dp.dacChannelOut` | The DAC output channel used for controlling the reward system. |

## `p.audio`

This structure contains all parameters related to auditory stimuli and feedback.

| Field | Description |
|---|---|
| `audsplfq` | The audio playback sampling rate in Hz for the DataPixx (e.g., 48000). |
| `Hitfq` | The frequency in Hz of the "high" tone, used to indicate a correct trial (hit). |
| `Missfq` | The frequency in Hz of the "low" tone, used to indicate an incorrect trial (miss/error). |
| `auddur` | The duration of the auditory tones in samples (e.g., 4800 samples, which is 100ms at 48kHz). |
| `lineOutLevel` | The audio level for the DataPixx line out, on a scale from 0 to 1. |
| `pcPlayback` | A boolean flag to determine if audio should be played from the host PC's soundcard (via Psychtoolbox) instead of the DataPixx. |

## `p.draw`

This structure defines parameters for visual elements that are drawn on the screen.

| Field | Description |
|---|---|
| `eyePosWidth` | The width in pixels of the eye position indicator on the experimenter's screen. |
| `fixPointWidth` | The line width in pixels of the fixation point. |
| `fixPointRadius` | The radius in pixels of the fixation point. |
| `fixWinPenPre` | The line width of the fixation window *before* a change. |
| `fixWinPenPost` | The line width of the fixation window *after* a change. |
| `fixWinPenDraw` | This variable is assigned the value of either `fixWinPenPre` or `fixWinPenPost` during the trial to dynamically change the fixation window's appearance. |
| `gridSpacing` | The spacing of the reference grid on the experimenter's display in degrees of visual angle. |
| `gridW` | The line width of the reference grid lines. |
| `joyRect` | The position and dimensions of the rectangle used to indicate joystick status on the experimenter's display. |
| `cursorW` | The width of the cursor in pixels. |
| `clutIdx` | This sub-structure contains indices for the Color Look-Up Table (CLUT). Each field name is a human-readable identifier for a color, and its value corresponds to a row in the CLUT. The format is `expColor_subColor`, where `exp` is the color on the experimenter's screen and `sub` is the color on the subject's screen. Example Fields: `expBlack_subBlack`, `expBg_subBg`, `expRed_subBg`, `expCyan_subCyan`, etc. |
| `color` | This sub-structure assigns specific CLUT indices from `p.draw.clutIdx` to different visual elements in the task. This is where the color of each component is defined for the current trial. |
| `color.background` | CLUT index for the screen's background color. |
| `color.cursor` | CLUT index for the cursor color. |
| `color.fix` | CLUT index for the fixation point color. |
| `color.fixWin` | CLUT index for the fixation window color. |
| `color.cueDots` | CLUT index for the cue dots color. |
| `color.foilDots` | CLUT index for the foil dots color. |
| `color.eyePos` | CLUT index for the eye position indicator. |
| `color.gridMajor` | CLUT index for the major grid lines on the experimenter screen. |
| `color.gridMinor` | CLUT index for the minor grid lines on the experimenter screen. |
| `color.cueRing` | CLUT index for the cue ring color. |
| `color.joyInd` | CLUT index for the joystick indicator color. |

## `p.state`

This structure defines the integer codes for the different states of the trial's state machine.

| State | Code | Description |
|---|---|---|
| `trialBegun` | 1 | The initial state of every trial. Used for setup before the ITI begins. |
| `waitForITI` | 2 | The inter-trial interval (ITI) state, a pause between trials. |
| `showCue` | 3 | The state where the fixation cue is displayed, prompting the subject to fixate. |
| `waitForFix` | 4 | The state where the system waits for the subject's gaze to enter the fixation window. |
| `holdFix` | 5 | The state where the subject must maintain fixation within the window for a specified duration. |
| `showOutcome` | 6 | The state where the token stimuli are displayed on the screen. |
| `cashInTokens` | 7 | The state where the token "cashing in" animation occurs, and rewards are delivered sequentially. |
| `fixBreak` | 11 | The trial is aborted because the subject looked away from the fixation point during the `holdFix` state. |
| `nonStart` | 12 | The trial is aborted because the subject failed to acquire fixation on the cue within the allotted time. |
| `success` | 21 | The trial was completed successfully, and all rewards were delivered. |

## `p.trVars`

This structure manages variables that can change on a trial-by-trial basis.

| Field | Description |
|---|---|
| `passJoy` | If `true`, simulates correct joystick trials for debugging. |
| `passEye` | If `true`, simulates correct eye position trials for debugging. |
| `blockNumber` | The current block number. |
| `repeat` | If `true`, the current trial's parameters will be repeated on the next trial. |
| `rwdJoyPR` | If `0`, reward is given for a joystick *press*. If `1`, reward is given for a joystick *release*. |
| `finish` | The number of trials after which the experiment will automatically stop. |
| `filesufix` | A suffix for the data file. |
| `joyVolt`, `eyeDegX`, `eyeDegY`, `eyePixX`, `eyePixY` | Live readings of joystick voltage and eye position in degrees and pixels. Initialized to 0. |
| `fixDegX`, `fixDegY` | The X and Y coordinates of the fixation point in degrees. |
| `fixDur` | The required duration for which the subject must hold fixation in seconds. |
| `fixAqDur` | The maximum time allowed to acquire fixation on the cue in seconds. |
| `fixWinWidthDeg`, `fixWinHeightDeg` | The width and height of the fixation window in degrees. |
| `rewardDurationMs` | The duration of a single juice pulse for a reward in milliseconds. |
| `juicePause` | The pause between sequential juice rewards when "cashing in" multiple tokens in seconds. |
| `outcomeDelay` | The delay after a successful fixation before the tokens are shown and "cashed in" in seconds. |
| `tokenI` | An index used to count through tokens during the reward delivery sequence. |
| `itiMean`, `itiMin`, `itiMax` | Parameters defining the duration of the inter-trial interval (ITI) in seconds. |
| `tokenBaseX`, `tokenBaseY` | The (X,Y) position of the first token in degrees. |
| `tokenSpacing` | The spacing between adjacent tokens in degrees. |
| `flickerFramesPerColor` | The number of screen refreshes each color is displayed for during the token flicker animation. |
| `currentState` | The current state of the state machine, initialized to `p.state.trialBegun`. |
| `exitWhileLoop` | A flag that, when set to `true`, terminates the `while` loop in the `_run.m` file, ending the current trial. |
| `fixPointRadPix`, `fixPointLinePix` | The radius and line weight of the fixation point in pixels. |
| `useCellsForDraw`, `wantEndFlicker`, `wantOnlinePlots` | Flags to control various drawing and plotting options. |
| `postFlip` | A structure used internally for precise timing, to log the exact time an event occurred after a screen flip. |
| `flipIdx` | An index that counts the number of screen flips in a trial. |

## `p.trData`

This structure contains all the data that is recorded during a single trial.

### Raw Data Streams
| Field | Description |
|---|---|
| `eyeX`, `eyeY` | Arrays containing the raw X and Y eye position data for the trial. |
| `eyeP` | Array containing the raw pupil diameter data. |
| `eyeT` | Array containing the timestamps for each eye data sample. |
| `joyV` | Array containing the raw joystick voltage readings. |
| `dInValues` | Array of values from the digital input channels on the DataPixx. |
| `dInTimes` | Array of timestamps for the digital input values. |
| `onlineEyeX`, `onlineEyeY` | Eye position data used for online plotting during the experiment. |

### Timing Information (`p.trData.timing`)
This sub-structure contains timestamps for all key events within the trial. Before each trial, all timing fields are initialized to -1. During the trial, they are populated with timestamps.

**All timestamps are measured in seconds relative to `p.trData.timing.trialStartPTB`**. This reference time is captured at the very beginning of the trial's `_run.m` function using `pds.getTimes`, which calls Psychtoolbox's `GetSecs`. Subsequent event times are calculated as `timeNow = GetSecs - p.trData.timing.trialStartPTB;`.

| Field | Description |
|---|---|
| `trialStartPTB` | The time the trial began, according to Psychtoolbox's `GetSecs`. This is the master reference time for all other timestamps in this structure. |
| `trialStartDP` | The time the trial began, according to the DataPixx clock. Can be used for synchronization. |
| `trialEnd` | The time the trial's `run` loop concluded. |
| `lastFrameTime` | The timestamp of the most recent screen flip, relative to `trialStartPTB`. |
| `flipTime` | An array containing the timestamp for every screen flip that occurred during the trial, relative to `trialStartPTB`. |
| `cueOn` | The time the fixation cue appeared on the screen. |
| `fixAq` | The time the subject's gaze first entered the fixation window. |
| `fixBreak` | The time the subject's gaze left the fixation window during the hold period. |
| `outcomeOn` | The time the token outcome was displayed. |
| `reward` | The time of reward delivery. This is a vector that stores the time of each reward pulse. |
| `stimOn` | The time the stimulus appeared (Note: this seems to be defined but not used in the `tokens_run.m` state machine, which uses `outcomeOn`). |
| `stimOff` | The time the stimulus was removed (Note: seems unused). |

## `p.stim`

This structure contains parameters related to the visual stimuli used in the task.

| Field | Description |
|---|---|
| `token.radius` | The radius of each token stimulus in degrees of visual angle. |
| `token.color` | The color of the token stimuli, defined as an [R, G, B] triplet. |
