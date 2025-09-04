function codes = initCodes
%   codes = pds.initCodes
%
% PAN-TASK function that initializes codes used to strobe events to the 
% ehpys recording system.
% These are the same codes that will identify events in the ephys file
% so this file is HOLY. 
% Once recording has been done, this file is the only way to reconstruct 
% the data so I'm not kidding-- H O L Y.
%
% For every new taks you code, it will likely use many event codes that are
% already present in this file. Enjoy them. For any new codes you might
% need, just add them to this file and verify that don't overlap with
% existing codes by running the verification cell at the bottom. 

% two kinds of strobes:
% Paired-strobe
% timing-strobe

%% task code
% Paired-strobe. Its pair is a unique task code that is set in the
% settings file, and takes the value of one of the unique task codes 
% defined in the cell "holy unique task codes", below.
codes.taskCode          = 32000;

%% holy unique task codes:
% Each task gets its own unique task code for easy identification. These
% are the values that are strobed after taskCode is strboed. 
codes.uniqueTaskCode_mcd        	= 32001;
codes.uniqueTaskCode_gSac    		= 32002;
codes.uniqueTaskCode_freeView   	= 32003;
codes.uniqueTaskCode_pFix       	= 32004;
codes.uniqueTaskCode_pFixLfp    	= 32005;
codes.uniqueTaskCode_pFixMotDir 	= 32006;
codes.uniqueTaskCode_mFlash     	= 32007;
codes.uniqueTaskCode_tod        	= 32008;
codes.uniqueTaskCode_scd        	= 32009;
codes.uniqueTaskCode_nfl        	= 32010;
codes.uniqueTaskCode_gSac_jph  		= 32011;
codes.uniqueTaskCode_gSac_contrast  = 32012;
codes.uniqueTaskCode_seansFirstTask = 32013;
codes.uniqueTaskCode_tokens         = 32014;
codes.uniqueTaskCode_gSac_4factors  = 32015;

%% unique codes that are internal to the 'classyStrobe' function class
% (see pds.classyStrobe.m for more details)
ss                  = pds.classyStrobe;
codes.isCell        = ss.internalStrobeCodes.isCell;
codes.cellLength    = ss.internalStrobeCodes.cellLength;


%% currently using fst codse...

% trial codes
codes.trialBegin        = 30001; % The very beginning of a trial 
codes.trialEnd          = 30009; % The very end of a trial 
codes.connectPLX        = 11001; % ???
codes.trialCount        = 11002; 
codes.blockNumber       = 11003; 
codes.trialInBlock      = 11004;
codes.setNumber         = 11005;
codes.state             = 11008;
codes.trialCode         = 11009;
codes.trialType         = 11010;
codes.fileSufix         = 11011; 
codes.taskType          = 11099;
codes.goodTrialCount    = 11100;
codes.goodtrialornot    = 21101; % Gongchen Added on 2019/12/30 I think it is better than goodTrialCount
                                 % Also I think these codes are better
                                 % start from '2~9' rather than 1, because
                                 % the code.time_1hhmm can sometimes contaminate the code  
%%  date & time
% these codes have a '1' before the date/time signifiers because a given
% date could lose its 0, e.g. the time_hhmm: 0932, would be sent as 932. By
% adding a '1' we get 10932, thus saving the 0. As long as user remembers
% to remove the first digit from the strobed values, we're all good.
codes.date_1yyyy      = 11102;
codes.date_1mmdd      = 11103;
codes.time_1hhmm      = 11104;

%%
codes.repeat20          = 11098; % 1 = 20 repeat trials during MemSac task.
codes.vissac            = 11097; % 1 = vis sac; 0 = memsac protocol
codes.inactivation      = 11095; % during inactivation
codes.useMotionStim     = 11094; % use motion stim for mapping


%% end of trial codes:
% code to represent a trial non strat 
codes.nonStart              = 22004;
codes.joyBreak              = 2005;
codes.fixBreak          = 3005;
codes.fixBreak2         = 3006; % this is if monkey breaks fixation whennot holding joystick in attn task

codes.saccToTargetOne	= 3007; % Used to identify which target monkey made a saccade to
codes.saccToTargetTwo	= 3008;

%% optical stimulation codes
codes.optoStimOn        = 17001;
codes.optoStimTrial     = 17002;
codes.optoStimSham      = 17003;

%% joystick codes
codes.joyPress              = 2001;
codes.joyRelease            = 2002;
codes.joyPressVoltDir       = 2010; 
codes.passJoy               = 2011;

%% fixation codes
codes.fixOn             = 3001;
codes.fixDim            = 3002;
codes.fixOff            = 3003;
codes.fixAq             = 3004;
codes.fixTheta          = 13001;
codes.fixRadius         = 13002;
codes.fixDimValue       = 13003;
codes.fixChangeTrial    = 13004;

%% saccade codes (used in gSac)
codes.saccadeOnset      = 2003;
codes.saccadeOffset     = 2004;
codes.blinkDuringSac    = 2007;

%% target codes (used in gSac)
codes.targetOn          = 4001;
codes.targetDim         = 4002;
codes.targetOff         = 4003;
codes.targetAq          = 4004;
codes.targetFixBreak    = 4005;
codes.targetReillum     = 4006; % target reillumination after a successful memory guided saccade
codes.targetTheta       = 14001;
codes.targetRadius      = 14002;

%% cue codes (used in mcd)
codes.cueOn             = 5001;
codes.cueOff            = 5003;
codes.stimLoc1Elev      = 15001;
codes.stimLoc1Ecc       = 15002;
codes.stimLoc2Elev      = 15003;
codes.stimLoc2Ecc       = 15004;

%% stimulus codes (used in mcd, pFix, etc.)

codes.stimOnDur                 = 5991;
codes.stimOffDur                = 5992;
codes.stimOn                    = 6002; % timing
codes.stimOff                   = 6003; % timing

codes.cueChange                 = 6004;
codes.foilChange                = 6005;
codes.noChange                  = 6006;
codes.isCueChangeTrial          = 6007;
codes.isFoilChangeTrial         = 6008;
codes.isNoChangeTrial           = 6009;
codes.cueMotionDelta            = 6010;
codes.foilMotionDelta           = 6011;
codes.cueStimIsOn               = 6012; % cued stimulus was shown in this trial
codes.foilStimIsOn              = 6013; % foil stimulus was shown in this trial
codes.isContrastChangeTrial     = 6014; % this trial had a contrast change
codes.hit                       = 6015; % this trial ended in a hit
codes.miss                      = 6016; % this trial ended in a miss
codes.foilFa                    = 6017; % this trial ended in a foil FA
codes.cr                        = 6018; % this trial ended in a CR
codes.fa                        = 6019; % this trial ended in a FA
codes.stimChange                = 6020;
codes.noChange                  = 6021;
codes.isStimChangeTrial         = 6022;
codes.isNoChangeTrial           = 6023;
codes.stimLoc1On                = 6024; % stimulus at location one was on in this trial
codes.stimLoc2On                = 6025; % stimulus at location one was on in this trial
codes.stimLoc3On                = 6026; % stimulus at location one was on in this trial
codes.stimLoc4On                = 6027; % stimulus at location one was on in this trial
codes.stimChangeTrial           = 16003;
codes.chgLoc                    = 16004;
codes.cueLoc                    = 16005;

% stimulus location & direction:
codes.stimLocRadius_x100  = 16001; % used to be named 'rfLocEcc'
codes.stimLocTheta_x10    = 16002; % used to be named 'rfLocTheta'
codes.stimMotDir          = 24000; % this is to send stim info; for eevnt codes for each dir see trialcodes.dirtun.m

% code for random number generation seeds
codes.stimSeed          = 16666;
codes.trialSeed         = 16667;

% code for orientations in orn tuning task
codes.orn               = 25000; % this is to send stim info; for eevnt codes for each orn see trialcodes.orntun.m

% reward code
codes.reward            = 8000;
codes.freeReward        = 8001;
codes.noFreeReward      = 8002;
codes.rewardDuration    = 18000;

% micro stim codes
codes.microStimOn       = 7001;

% audio codes
codes.audioFBKon        = 9000;
codes.lowTone           = 9001;
codes.noiseTone         = 9002;
codes.highTone          = 9003;

% image codes (used in freeview)
codes.imageId           = 6660;     % id of image
codes.imageOn           = 6661;     % time of image onset
codes.imageOff          = 6662;     % time of image offset
codes.freeViewDur       = 6663;     % duration of free image viewing

% tokens task codes
codes.CUE_ON = 5; % [cite: 1]
codes.REWARD_GIVEN = 7; % [cite: 1]
codes.TRIAL_END = 6; % [cite: 1]
codes.REWARD_AMOUNT_BASE = 100; % Base for strobing reward amount [cite: 1]
codes.OUTCOME_DIST_BASE = 90; % Base for strobing outcome distribution type [cite: 1]
codes.rwdAmt = 101;

% gSac_4factors codes
codes.halfBlock     = 16010;
codes.stimType      = 16011;
codes.salience      = 16012;
codes.targetColor   = 16013;
codes.targetLocIdx  = 16014;

%% validation

% making sure that every code is listed once (i.e. unique)
values = [];
flds = fieldnames(codes);
for iF = 1:numel(flds)
    if ~isstruct(codes.(flds{iF})) % so that we don't go over the uniqueTaskCode struct
        values(iF) = codes.(flds{iF}); %#ok<AGROW>
    end
end
% find duplicate values:
uValues         = unique(values);
ptr             = hist(values, uValues)>1;
duplicateCodes  = uValues(ptr);

if isempty(duplicateCodes)
    % all good! no duplicate values!
else
    disp('you have duplicate code numbers!!!')
    disp(duplicateCodes)
    error('YOU MUST FIX IT')
    keyboard
end

