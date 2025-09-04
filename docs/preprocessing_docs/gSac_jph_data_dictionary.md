# gSac_jph Saved Data Dictionary

This document provides a comprehensive description of the saved data structures for the 'gSac_jph' task. The overall 'p' structure is saved once at the beginning of a session, and this document explains all its subfields. The `p.trVars` and `p.trData` structures are saved on each trial.

## `p.init`

The `p.init` structure contains variables that are initialized once at the beginning of the experiment.

| Field | Description |
|---|---|
| `pcName` | The name of the computer running the experiment. |
| `rigConfigFile` | The path to the rig configuration file. This file contains subject/rig-specific details. |
| `taskName` | The name of the task (e.g., 'gSac_jph'). |
| `pldapsFolder` | The path to the PLDAPS folder. |
| `protocol_title` | The title of the experimental protocol, used for display purposes. |
| `date` | The date of the experiment in 'yyyymmdd' format. |
| `time` | The time the experiment started in 'HHMM' format. |
| `outputFolder` | The folder where the output files are saved. |
| `sessionId` | A unique identifier for the session, combining date, time, and task name. |
| `sessionFolder` | The folder for the current session's data. |
| `taskFiles` | A structure containing the names of the task's M-files: `init`, `next`, `run`, and `finish`. |
| `useDataPixxBool` | A boolean indicating whether the DataPixx/ViewPixx is being used. |
| `taskActions` | A cell array of strings with the names of M-files for user-defined actions accessible from the GUI. |
| `exptType` | A string indicating which version of the experiment is being run (e.g., 'all_locs', 'two_locs'). |
| `trDataInitList` | A cell array that lists trial data variables and their initial values. These are reset at the start of each trial. This list is defined in the `_settings.m` file and used by `initTrData.m`. |
| `nTrDataListRows`| The number of rows in `trDataInitList`. |
| `strobeList` | A cell array defining variables to be strobed to the ephys system at the end of each trial. |
| `codes` | A structure containing numeric codes for various trial events (e.g., `trialBegin`, `fixOn`). These are used for strobing. |
| `strb` | An instance of the `pds.classyStrobe` class, used to manage strobing of event codes. |

## `p.rig`

The `p.rig` structure contains rig-specific information, such as screen parameters, joystick calibration, and hardware settings. Many of these values are loaded from a rig-specific configuration file (e.g., `rigConfig_rig1.m`).

| Field | Description |
|---|---|
| `screen_number` | The screen number on which the experiment is displayed (e.g., 1 for an external monitor). |
| `refreshRate` | The refresh rate of the screen in Hz. |
| `frameDuration` | The duration of a single frame in seconds (1 / `refreshRate`). |
| `joyThreshPress` | The joystick voltage threshold for detecting a press. |
| `joyThreshRelease` | The joystick voltage threshold for detecting a release. |
| `magicNumber` | A small time offset (in seconds) used to ensure screen flips are not missed. |
| `joyVoltageMax` | The maximum voltage of the joystick. |
| `guiStatVals` | A cell array of strings listing status variables to be displayed in the GUI. |
| `guiVars` | A cell array of strings listing trial variables to be displayed and controlled from the GUI. |
| `viewdist` | The viewing distance from the subject's eyes to the screen in millimeters. |
| `screenhpix` | The height of the screen in pixels. |
| `screenh` | The height of the screen in millimeters. |
| `deg2PixConstant` | A constant used to convert degrees of visual angle to pixels. |
| `screen_name` | A name for the screen (e.g., 'viewpixx'). |
| `baseReward` | The base duration for reward delivery in milliseconds. |
| `joyVoltageMin` | The minimum voltage of the joystick. |
| `joyVoltageRest` | The resting voltage of the joystick. |
| `joyThresh` | A structure containing detailed joystick thresholds for high and low voltage presses and releases. |
| `dp` | A structure containing settings for the DataPixx/ViewPixx, such as ADC/DAC rates and buffer addresses. |
| `connectToOmniplex` | A boolean indicating whether to connect to an Omniplex system. |

## `p.audio`

The `p.audio` structure contains parameters for auditory feedback.

| Field | Description |
|---|---|
| `audsplfq` | The audio sampling frequency in Hz. |
| `Hitfq` | The frequency of the tone played on a correct trial in Hz. |
| `Missfq` | The frequency of the tone played on an incorrect trial in Hz. |
| `auddur` | The duration of the audio tone in samples. |

## `p.draw`

The `p.draw` structure contains parameters related to drawing visual stimuli on the screen.

| Field | Description |
|---|---|
| `clutIdx` | A structure that maps human-readable color names to indices in the Color Look-Up Table (CLUT). This provides a way to define colors centrally and refer to them by name. |
| `color` | A structure that holds the CLUT indices for different visual elements in the task (e.g., `p.draw.color.background`, `p.draw.color.fix`). The values are set from `p.draw.clutIdx`. |
| `fixPointWidth` | The width of the fixation point indicator line in pixels. |
| `fixPointRadius`| The "radius" of the fixation point in pixels. |
| `fixWinPenThin` | The pen width for drawing the fixation window before the subject is required to respond. |
| `fixWinPenThick`| The pen width for drawing the fixation window after the 'go' signal. |
| `fixWinPenDraw` | The current pen width for the fixation window, assigned either the thin or thick value during the trial. |
| `targWinPenThin`| The pen width for drawing the target window. |
| `targWinPenThick`| The pen width for drawing the target window. |
| `targWinPenDraw`| The current pen width for the target window. |
| `eyePosWidth` | The width of the eye position indicator in pixels. |
| `gridSpacing` | The spacing of the grid on the experimenter's display in degrees of visual angle. |
| `gridW` | The line width of the grid on the experimenter's display. |
| `joyRect` | A 4-element vector defining the rectangle for the joystick indicator on the experimenter's display. |
| `cursorW` | The width of the mouse cursor in pixels. |

## `p.state`

The `p.state` structure defines the different states of the trial's state machine. Each field corresponds to a state and holds a unique integer identifier.

| Field | Description |
|---|---|
| `trialBegun` | The initial state at the beginning of a trial. |
| `waitForJoy` | The state where the system is waiting for the subject to press and hold the joystick. |
| `showFix` | The state where the fixation point is shown, and the system waits for the subject to acquire fixation. |
| `dontMove` | The state where the subject must maintain fixation and hold the joystick. |
| `makeSaccade` | The state after the 'go' signal (fixation offset), where the subject is expected to make a saccade to the target. |
| `checkLanding` | The state where the system checks if the saccade has landed within the target window. |
| `holdTarg` | The state where the subject must maintain fixation on the target after a successful saccade. |
| `sacComplete` | A successful trial completion state, leading to reward. |
| `fixBreak` | An error state entered if the subject breaks fixation. |
| `joyBreak` | An error state entered if the subject releases the joystick prematurely. |
| `nonStart` | An error state entered if the subject fails to start the trial (e.g., by not pressing the joystick). |
| `failedToHoldTarg` | An error state if the subject fails to maintain fixation on the target. |

## `p.status`

The `p.status` structure contains variables that track the overall status and performance during the experiment.

| Field | Description |
|---|---|
| `iTrial` | The current trial number. |
| `iGoodTrial` | The number of successfully completed trials. |
| `iGoodVis` | The number of successful visually-guided saccade trials. |
| `iGoodMem` | The number of successful memory-guided saccade trials. |
| `pGoodVis` | The proportion of successful visually-guided saccade trials. |
| `pGoodMem` | The proportion of successful memory-guided saccade trials. |
| `iTarget` | An index or identifier for the current target. |
| `rippleOnline` | A flag indicating the status of the Ripple neurophysiology recording system. |

## `p.trVars`

The `p.trVars` structure holds all the variables that can change on a trial-by-trial basis. It inherits its base values from `p.trVarsInit` (defined in the settings file) and `p.trVarsGuiComm` (updated by the GUI), but many values, especially for timing and target location, are set dynamically in the `_next.m` file for each trial. This entire structure is saved for every trial.

### General & Control Variables

| Field | Description |
|---|---|
| `passJoy` | If set to 1, simulates correct joystick behavior. Useful for debugging. |
| `passEye` | If set to 1, simulates correct eye movement behavior. Useful for debugging. |
| `connectPLX` | A flag indicating whether to connect to the Plexon system. |
| `joyPressVoltDirection`| Defines the direction of joystick press (-1 for down, 1 for up). |
| `blockNumber` | The current block number. |
| `repeat` | A flag that, if true, indicates the trial should be repeated. |
| `rwdJoyPR` | Determines reward delivery: 0 for joystick press, 1 for release. |
| `wantEndFlicker` | If true, the screen flickers at the end of the trial if the joystick is still held. |
| `finish` | The maximum number of trials for the session. |
| `filesufix` | A numerical suffix for the saved data file. |
| `setTargLocViaMouse` | If true, the target location for the next trial is set by a mouse click. |
| `setTargLocViaGui` | If true, the target location can be set manually from the GUI. |
| `setTargLocViaTrialArray` | If true, target locations are drawn from a predefined list (`p.trVars.trialsArray`). |
| `connectRipple` | A boolean indicating whether to connect to the Ripple neurophysiology system. |
| `rippleChanSelect` | Selects a specific Ripple channel to use for online analysis. |
| `useOnlineSort` | If true, uses spike times sorted online in Trellis. |
| `wantOnlinePlots` | If true, online performance plots are generated and updated. |

### Stimulus, Geometry & Position Variables

| Field | Description |
|---|---|
| `propVis` | The proportion of trials that are visually-guided (as opposed to memory-guided). |
| `isVisSac` | A boolean that is true if the current trial is visually-guided. |
| `fixDegX`, `fixDegY` | The X and Y coordinates of the fixation point in degrees of visual angle. |
| `targDegX`, `targDegY` | The X and Y coordinates of the target in degrees of visual angle. |
| `targTheta`, `targRadius` | The polar coordinates (angle and eccentricity) of the target, calculated from `targDegX` and `targDegY`. |
| `targRadius_x100` | The target radius (eccentricity) multiplied by 100 for strobing. |
| `targTheta_x10` | The target angle (in degrees) multiplied by 10 for strobing. |
| `minTargAmp`, `maxTargAmp` | The minimum and maximum allowed target amplitude (eccentricity). |
| `staticTargAmp` | A fixed target amplitude, if used. |
| `maxHorzTargAmp`, `maxVertTargAmp` | Maximum horizontal and vertical target amplitudes for 'rectangular annulus' target selection. |
| `fixWinWidthDeg`, `fixWinHeightDeg` | The width and height of the fixation window in degrees. |
| `targWinWidthDeg`, `targWinHeightDeg` | The width and height of the target window in degrees. |
| `targWidth` | The line width of the target rectangle in pixels. |
| `targRadius` | This variable is overloaded. It is initialized as the target size (radius in pixels) in `trVarsInit`, but then overwritten in `nextParams` with the target's eccentricity in degrees. The drawing code uses this eccentricity value as the size in pixels. |
| `stim` | A copy of the `p.stim` structure for the current trial. |
| `joyVolt` | The current voltage reading from the joystick. |
| `eyeDegX`, `eyeDegY` | The current X and Y eye position in degrees. |
| `eyePixX`, `eyePixY` | The current X and Y eye position in pixels. |
| `mouseCursorX`, `mouseCursorY`| The current X and Y coordinates of the mouse cursor. |

### Timing & Duration Variables

| Field | Description |
|---|---|
| `rewardDurationMs` | The duration of the reward delivery in milliseconds. |
| `rewardDelay` | The delay between a correct response and the delivery of the reward. |
| `timeoutAfterFa` | The timeout duration after a false alarm. |
| `joyWaitDur` | The maximum time to wait for the joystick to be pressed at the start of a trial. |
| `fixWaitDur` | The maximum time to wait for the subject to acquire fixation after pressing the joystick. |
| `freeDur` | A duration at the beginning of the trial with no requirements. |
| `trialMax` | The maximum possible duration of a single trial. |
| `joyReleaseWaitDur`| The time to wait after the trial ends before starting the end-of-trial flicker. |
| `postRewardDuration` | The duration the trial continues after reward delivery to record neural activity. |
| `targetFlashDuration` | For memory-guided trials, the duration the target is visible before disappearing. |
| `targHoldDurationMin`, `targHoldDurationMax` | The minimum and maximum required duration for holding fixation on the target. |
| `targHoldDuration` | The specific target hold duration for the current trial, randomized between the min and max values. |
| `maxSacDurationToAccept` | The maximum allowed duration for a saccade to be considered valid. |
| `goLatencyMin`, `goLatencyMax` | The minimum and maximum allowed reaction time for initiating a saccade after the 'go' signal. |
| `targOnsetMin`, `targOnsetMax` | The minimum and maximum time from fixation acquisition to target onset. |
| `timeTargOnset` | The specific time from fixation acquisition to target onset for the current trial. |
| `timeTargOffset` | The specific time from fixation acquisition to target offset for the current trial. |
| `goTimePostTargMin`, `goTimePostTargMax` | The minimum and maximum time from target onset to the 'go' signal (fixation offset). |
| `timeFixOffset` | The specific time from fixation acquisition to the 'go' signal for the current trial. |
| `maxFixWait` | The maximum time allowed to acquire fixation. |
| `targOnSacOnly` | If true, the target in memory trials only reappears contingent on a saccade. |
| `rwdTime` | The time of reward delivery, initialized to -1. |
| `targTrainingDelay` | A delay for target re-illumination during training stages. |
| `timeoutdur` | The duration of the timeout period after an error. |

### Internal State & Loop Variables

| Field | Description |
|---|---|
| `currentState` | A numeric value representing the current state of the trial's state machine. |
| `exitWhileLoop` | A boolean flag that, when set to true, terminates the main `while` loop of the trial. |
| `targetIsOn` | A boolean flag indicating if the target is currently being displayed. |
| `postMemSacTargOn` | A flag to turn the target on after a memory-guided saccade is completed. |
| `whileLoopIdx` | A counter for the number of iterations of the main `while` loop in the `_run` function. |
| `flipIdx` | A counter for the number of screen flips in the trial. |
| `stimFrameIdx` | A counter for stimulus frames. |
| `eyeVelFiltTaps` | The number of samples used for the online eye velocity filter. |
| `eyeVelThresh` | The velocity threshold (in deg/s) for online saccade detection. |
| `useVelThresh` | A boolean to enable or disable the use of the velocity threshold for online saccade detection. |
| `eyeVelThreshOffline` | The velocity threshold used for offline saccade analysis. |

## `p.trData`

The `p.trData` structure is where all the data collected during a single trial is stored. It is initialized as an empty structure at the beginning of each trial and populated throughout the trial's execution. This entire structure is saved for every trial.

### Raw Data Buffers

| Field | Description |
|---|---|
| `eyeX`, `eyeY` | Buffers containing the raw X and Y eye position data from the eye tracker. |
| `eyeP` | Buffer containing the raw pupil size data. |
| `eyeT` | Buffer containing the timestamps for the eye data samples. |
| `joyV` | Buffer containing the raw joystick voltage readings. |
| `dInValues`, `dInTimes` | Buffers for digital input values and their corresponding timestamps. |
| `onlineGaze` | A matrix containing online gaze calculations: [X (deg), Y (deg), time (s), velocity (deg/s)]. |
| `strobed` | A list of all the numerical codes that were strobed to the ephys system during the trial. |
| `spikeTimes`, `spikeClusters` | Data from the Ripple system: timestamps and cluster IDs for detected spikes. |
| `eventTimes`, `eventValues` | Data from the Ripple system: timestamps and values for recorded events. |

### Calculated Saccade & Trial Parameters

| Field | Description |
|---|---|
| `preSacXY` | The average X,Y eye position in the 5ms before saccade onset. |
| `postSacXY` | The average X,Y eye position in the 5ms after saccade offset. |
| `peakVel` | The peak velocity of the primary saccade. |
| `SRT` | Saccade Reaction Time, calculated as the time from the 'go' signal to saccade onset. |
| `trialEndState` | The numerical identifier of the state in which the trial ended. |
| `trialRepeatFlag` | A boolean flag that is true if the trial was an error and should be repeated. |

### Timing Information

The `p.trData.timing` substructure contains timestamps for all critical events within the trial. Before each trial, all timing fields are initialized to -1. During the trial, they are populated with timestamps.

**All timestamps are measured in seconds relative to `p.trData.timing.trialStartPTB`**. This reference time is captured at the very beginning of the trial's `_run.m` function using `pds.getTimes`, which calls Psychtoolbox's `GetSecs`. Subsequent event times are calculated as `timeNow = GetSecs - p.trData.timing.trialStartPTB;`.

| Field | Description |
|---|---|
| `trialStartPTB` | The trial start time as recorded by Psychtoolbox (`GetSecs`). This is the master reference time for all other timestamps in this structure. |
| `trialStartDP` | The trial start time as recorded by the DataPixx. Can be used for synchronization. |
| `trialBegin` | Timestamp of when the trial logic began (state 1: `trialBegun`). |
| `joyPress` | Timestamp of when the joystick was pressed to initiate the trial. |
| `trialEnd` | Timestamp of when the trial logic ended (when `exitWhileLoop` becomes true). |
| `lastFrameTime` | The timestamp of the last screen flip, relative to `trialStartPTB`. |
| `flipTime` | A vector containing the timestamps of every screen flip during the trial, relative to `trialStartPTB`. |
| `fixOn` | Timestamp of when the fixation point appeared. |
| `fixAq` | Timestamp of when fixation was acquired. |
| `fixOff` | Timestamp of when the fixation point disappeared (the 'go' signal). |
| `targetOn` | Timestamp of when the saccade target appeared. |
| `targetOff` | Timestamp of when the saccade target disappeared (in memory-guided trials). |
| `targetReillum` | Timestamp of when the target reappeared (in memory-guided trials). |
| `targetAq` | Timestamp of when the target was acquired (i.e., when the saccade landed in the target window). |
| `saccadeOnset` | Timestamp of when the saccade was initiated (detected by eye leaving fixation window or by velocity threshold). |
| `saccadeOffset` | Timestamp of when the saccade ended (detected by eye entering target window). |
| `brokeFix` | Timestamp of when a fixation break occurred. |
| `joyRelease` | Timestamp of when the joystick was released prematurely. |
| `reward` | Timestamp of when the reward was delivered. |
| `tone` | Timestamp of when the auditory feedback was delivered. |
| `frameNow` | The current frame number within the trial, calculated as `fix((GetSecs - trialStartPTB) * refreshRate)`. |

## `p.stim`

The `p.stim` structure holds parameters that define the stimuli, primarily the set of possible target locations. A full copy of this structure is saved within `p.trVars` on every trial.

| Field | Description |
|---|---|
| `targLocationPreset` | A string specifying the method for generating target locations ('grid', 'ring', or 'nRing'). |
| `grid*` | If `targLocationPreset` is 'grid', these fields (`gridMinX`, `gridMaxX`, `gridBinSizeX`, etc.) define the boundaries and spacing of the target grid. |
| `ring*` | If `targLocationPreset` is 'ring' or 'nRing', these fields (`ringRadius`, `ringTargNumber`, `ringBaseAngle`) define the radius, number of targets, and starting angle for the target ring(s). |
| `dotWidth` | The width of the target dot stimulus in pixels. |
| `targList` | An N-by-2 matrix containing the [X, Y] coordinates (in degrees) of all possible target locations, generated based on the `targLocationPreset`. |
