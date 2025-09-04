function comps = local_defineComparisons()
% defines trial condition comparisons for the 4-factor gSac task and the attention task
% 4f gSac factors: stimulus type (face/nonface/bullseye), salience (hi/lo; bullseye only),
%                  probability (hi/lo), reward value via hemifield (hi/lo), location (t1/t2)

%% gSac (4-factor) comparisons — single-factor only
i = 0; comps.gSac = struct([]);

% target location
i=i+1; comps.gSac(i).name = 'targetLoc';
comps.gSac(i).condition1 = 'isT1Trial';
comps.gSac(i).condition2 = 'isT2Trial';
comps.gSac(i).description = 'target location (t1 vs t2)';

% stimulus type: face vs non-face
i=i+1; comps.gSac(i).name = 'stim_face_vs_nonface';
comps.gSac(i).condition1 = 'isFaceTrial';
comps.gSac(i).condition2 = 'isNonFaceTrial';
comps.gSac(i).description = 'stimulus type (face vs non-face)';

% stimulus type: face vs bullseye
i=i+1; comps.gSac(i).name = 'stim_face_vs_bullseye';
comps.gSac(i).condition1 = 'isFaceTrial';
comps.gSac(i).condition2 = 'isBullseyeTrial';
comps.gSac(i).description = 'stimulus type (face vs bullseye)';

% probability: high vs low
i=i+1; comps.gSac(i).name = 'probability';
comps.gSac(i).condition1 = 'isHighProbTrial';
comps.gSac(i).condition2 = 'isLowProbTrial';
comps.gSac(i).description = 'location probability (high vs low)';

% salience (bullseye only): high vs low
i=i+1; comps.gSac(i).name = 'salience_bullseye_only';
comps.gSac(i).condition1 = 'isHighSaliencyTrial';
comps.gSac(i).condition2 = 'isLowSaliencyTrial';
comps.gSac(i).description = 'salience (high vs low; bullseye only)';

% reward via hemifield: high vs low
i=i+1; comps.gSac(i).name = 'reward';
comps.gSac(i).condition1 = 'isHighRewardTrial';
comps.gSac(i).condition2 = 'isLowRewardTrial';
comps.gSac(i).description = 'reward value via hemifield (high vs low)';

%% attention comparisons (unchanged)
comps.attn = struct([]);
comps.attn(1).name = 'cueLoc';
comps.attn(1).condition1 = 'isCueInTrial';
comps.attn(1).condition2 = 'isCueOutTrial';
comps.attn(1).description = 'cue location (in vs out rf)';

comps.attn(2).name = 'cueChgLoc';
comps.attn(2).condition1 = 'isCueInChgInTrial';
comps.attn(2).condition2 = 'isCueOutChgOutTrial';
comps.attn(2).description = 'cue change location (in vs out rf, change at cued loc)';

comps.attn(3).name = 'cueChgInHitVsMiss';
comps.attn(3).condition1 = 'isCueInChgInHitTrial';
comps.attn(3).condition2 = 'isCueInChgInMissTrial';
comps.attn(3).description = 'cue change in rf (hit vs miss)';
end
