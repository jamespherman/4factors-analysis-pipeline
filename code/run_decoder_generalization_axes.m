function results = run_decoder_generalization_axes()
% RUN_DECODER_GENERALIZATION_AXES  SVM generalization + two-null axis test.
%
%   For each SC session, decode a target feature from target-onset
%   population activity at location 1, and ask whether the decoder (and
%   its coding axis) built under one reward/probability condition
%   transfers to the other.
%
%   Two decoded features (separate analyses):
%       salience  : high vs low salience,  BULLSEYE trials @ loc 1
%       identity  : face vs non-face,      IMAGE    trials @ loc 1
%
%   Two generalization dimensions (binary):
%       reward      : high (>200) vs low (<200) reward   [pooled over block]
%       probability : block 1 (high prob) vs block 2 (low prob) @ loc 1
%
%   Two population variants per session:
%       classic     : classic SC neurons only
%       ci          : classic + interneuron pooled
%
%   Activity: mean firing rate per neuron in the 50-250 ms target-onset
%   window (condition_defs...window_roc_plan.epoch_windows.visual.SC).
%   All successful location-1 trials are used (no stability filter: the
%   per-neuron stable-period intersection is empty for these populations).
%
%   GENERALIZATION: train a linear SVM on one condition, test on the
%   other (both directions), vs within-condition 5-fold CV, vs a
%   label-shuffle chance null.
%
%   TWO-NULL AXIS TEST (after Herman-lab color-contrast method):
%       observed : bootstrap-resample trials in each condition, train an
%                  SVM in each, cosine similarity of the weight vectors.
%       relabel  : shuffle the decoded-feature labels within each
%                  condition -> null for "is there any real axis".
%       regroup  : pool trials by feature level across the two
%                  conditions, re-split into two pseudo-conditions ->
%                  null for "the two conditions share ONE common axis".
%   Per-session verdict:
%       Inconclusive     : observed not above the relabel null.
%       Indistinguishable: above relabel AND overlaps regroup null.
%       Different        : above relabel AND below regroup null.
%
%   NOTE (probability confound): at location 1, probability = block, and
%   blocks are temporally separated, so the probability generalization is
%   confounded with recording drift / time-on-task. Interpret the
%   probability results with that in mind; the reward split is
%   interleaved within blocks and is not subject to this confound.

%% Setup
rng(42);
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
[~, cdefs] = define_task_conditions();
win = cdefs.window_roc_plan.epoch_windows.visual.SC;   % [0.05 0.25]

codes = initCodes();
man = readtable(fullfile(project_root, 'config', 'session_manifest.csv'));
sc = man(strcmp(man.brain_area, 'SC'), :);

%% Parameters
P = struct();
P.window       = win;
P.min_per_class = 8;    % min trials per class per condition
P.min_neurons   = 5;    % min neurons to build an axis
P.n_cv          = 5;    % folds for within-condition CV
P.n_perm        = 200;  % label-shuffle permutations (chance)
P.n_boot        = 200;  % bootstrap reps (observed cosine)
P.n_relabel     = 400;  % relabel-null reps
P.n_regroup     = 400;  % regroup-null reps
P.alpha         = 0.05;

variants = {'classic', 'ci'};
factors  = struct(...
    'name', {'salience', 'identity'}, ...
    'trialset', {'bullseye', 'image'}, ...
    'levelA', {'high', 'face'}, ...
    'levelB', {'low', 'nonface'});
genDims  = {'reward', 'probability'};

results = struct();
results.params = P;
results.sessions = {};

%% Per-session loop
for k = 1:height(sc)
    sid = sc.unique_id{k};
    f = fullfile(project_root, 'data', 'processed', sid, ...
        [sid '_session_data.mat']);
    if ~exist(f, 'file'); continue; end
    S = load(f, 'session_data'); sd = S.session_data;

    si = sd.analysis.screening_info;
    nclass = {si.neuron_class};

    % Window-averaged rates for ALL neurons: [nAllNeurons x nTrials]
    ev = sd.analysis.core_data.spikes.targetOn;
    tv = ev.time_vector;
    bmask = (tv >= P.window(1)) & (tv <= P.window(2));
    wrate = mean(ev.rates(:, :, bmask), 3, 'omitnan');  % [nNeurons x nTrials]

    % Trial condition vectors aligned to core_data columns
    tidx = sd.analysis.core_data.trial_indices;
    ti = sd.trialInfo;
    loc1 = (ti.targetLocIdx(tidx) == 1);
    hr = ti.rewardDuration(tidx) > 200;  lr = ti.rewardDuration(tidx) < 200;
    b1 = ti.blockNumber(tidx) == 1;      b2 = ti.blockNumber(tidx) == 2;
    bull = ti.stimType(tidx) > 2;        img = ismember(ti.stimType(tidx), [1 2]);
    lvl = struct();
    lvl.salience = struct('A', ti.salience(tidx) == 1, 'B', ti.salience(tidx) == 2, ...
                          'set', bull);
    lvl.identity = struct('A', ti.stimType(tidx) == 1, 'B', ti.stimType(tidx) == 2, ...
                          'set', img);
    gsel = struct('reward', {{hr, lr}}, 'probability', {{b1, b2}});

    sres = struct('session_id', sid);

    for iv = 1:numel(variants)
        variant = variants{iv};
        if strcmp(variant, 'classic')
            rows = find(strcmp(nclass, 'classic'));
        else
            rows = find(ismember(nclass, {'classic', 'interneuron'}));
        end
        Xall = wrate(rows, :)';   % [nTrials x nNeurons]

        for ifac = 1:numel(factors)
            fac = factors(ifac);
            base = loc1 & lvl.(fac.name).set;
            yA = lvl.(fac.name).A;   % +1 level
            yB = lvl.(fac.name).B;   % -1 level

            for ig = 1:numel(genDims)
                gd = genDims{ig};
                g1 = gsel.(gd){1};   g2 = gsel.(gd){2};

                key = sprintf('%s_%s_%s', variant, fac.name, gd);

                r = analyze_one(Xall, base, yA, yB, g1, g2, P);
                r.n_neurons = numel(rows);
                sres.(key) = r;
            end
        end
    end
    results.sessions{end+1} = sres;
    fprintf('DECODE %-24s done (classic=%d, ci=%d)\n', sid, ...
        sum(strcmp(nclass,'classic')), ...
        sum(ismember(nclass,{'classic','interneuron'})));
end

%% Save
outfile = fullfile(project_root, 'data', 'processed', ...
    'decoder_generalization_axes.mat');
save(outfile, 'results', '-v7.3');
fprintf('DECODE saved: %s\n', outfile);

end

%% ========================================================================
function r = analyze_one(Xall, base, yA, yB, g1, g2, P)
% One (variant x factor x genDim) cell: generalization + two-null axis.

r = struct('ok', false, 'reason', '', 'verdict', 'skipped');

if size(Xall, 2) < P.min_neurons
    r.reason = 'too_few_neurons'; return;
end

% Condition 1 / 2 trial index sets (within base), with labels
c1A = find(base & g1 & yA);  c1B = find(base & g1 & yB);
c2A = find(base & g2 & yA);  c2B = find(base & g2 & yB);

n = [numel(c1A) numel(c1B) numel(c2A) numel(c2B)];
r.n_trials = n;
if min(n) < P.min_per_class
    r.reason = 'too_few_trials'; return;
end

% Drop trials with any NaN neuron (missing event) up front
keep = @(idx) idx(~any(isnan(Xall(idx, :)), 2));
c1A = keep(c1A); c1B = keep(c1B); c2A = keep(c2A); c2B = keep(c2B);
if min([numel(c1A) numel(c1B) numel(c2A) numel(c2B)]) < P.min_per_class
    r.reason = 'too_few_trials_after_nan'; return;
end

X1 = Xall([c1A; c1B], :);  y1 = [ones(numel(c1A),1); -ones(numel(c1B),1)];
X2 = Xall([c2A; c2B], :);  y2 = [ones(numel(c2A),1); -ones(numel(c2B),1)];

svm = @(X, y) fitcsvm(X, y, 'KernelFunction', 'linear', 'Standardize', true);

%% --- Generalization ---
% cross-condition (both directions)
m1 = svm(X1, y1);  acc_1to2 = mean(predict(m1, X2) == y2);
m2 = svm(X2, y2);  acc_2to1 = mean(predict(m2, X1) == y1);
% within-condition CV
acc_cv1 = cv_acc(X1, y1, P.n_cv);
acc_cv2 = cv_acc(X2, y2, P.n_cv);
% label-shuffle chance null for cross accuracy (both directions pooled)
perm_acc = zeros(P.n_perm, 1);
for p = 1:P.n_perm
    ys1 = y1(randperm(numel(y1)));
    ys2 = y2(randperm(numel(y2)));
    mp1 = svm(X1, ys1); a12 = mean(predict(mp1, X2) == ys2);
    mp2 = svm(X2, ys2); a21 = mean(predict(mp2, X1) == ys1);
    perm_acc(p) = (a12 + a21) / 2;
end
obs_cross = (acc_1to2 + acc_2to1) / 2;
r.gen = struct('acc_1to2', acc_1to2, 'acc_2to1', acc_2to1, ...
    'acc_cross_mean', obs_cross, 'acc_cv1', acc_cv1, 'acc_cv2', acc_cv2, ...
    'perm_mean', mean(perm_acc), 'perm_p', ...
    (sum(perm_acc >= obs_cross) + 1) / (P.n_perm + 1));

%% --- Two-null axis test ---
obsCos = axis_cos_boot(Xall, c1A, c1B, c2A, c2B, P.n_boot, svm, false);
signFlip = median(obsCos) < 0;
if signFlip
    obsCos = axis_cos_boot(Xall, c1A, c1B, c2A, c2B, P.n_boot, svm, true);
end
relCos = axis_cos_relabel(Xall, c1A, c1B, c2A, c2B, P.n_relabel, svm, signFlip);
regCos = axis_cos_regroup(Xall, c1A, c1B, c2A, c2B, P.n_regroup, svm);

medObs = median(obsCos);
p_vs_relabel = mean(relCos >= medObs);           % small => above relabel
regroup_lo   = prctile(regCos, 100 * P.alpha);   % lower tail of regroup
above_relabel = p_vs_relabel < P.alpha;
overlaps_regroup = medObs >= regroup_lo;

if ~above_relabel
    verdict = 'Inconclusive';
elseif overlaps_regroup
    verdict = 'Indistinguishable';
else
    verdict = 'Different';
end

r.axis = struct('obsCos', obsCos, 'relabelCos', relCos, 'regroupCos', regCos, ...
    'median_obs', medObs, 'p_vs_relabel', p_vs_relabel, ...
    'regroup_lo', regroup_lo, 'signFlip', signFlip);
r.verdict = verdict;
r.ok = true;

end

%% ------------------------------------------------------------------------
function a = cv_acc(X, y, nfold)
% Stratified k-fold CV accuracy for a linear standardized SVM.
try
    c = cvpartition(y, 'KFold', min(nfold, min(sum(y==1), sum(y==-1))));
catch
    a = NaN; return;
end
pred = nan(size(y));
for i = 1:c.NumTestSets
    tr = training(c, i); te = test(c, i);
    m = fitcsvm(X(tr,:), y(tr), 'KernelFunction','linear','Standardize',true);
    pred(te) = predict(m, X(te,:));
end
a = mean(pred == y);
end

%% ------------------------------------------------------------------------
function cosv = axis_cos_boot(Xall, c1A, c1B, c2A, c2B, nrep, svm, flip)
% Bootstrap observed cosine: resample within each class of each condition.
cosv = zeros(nrep, 1);
for b = 1:nrep
    i1A = c1A(randi(numel(c1A), numel(c1A), 1));
    i1B = c1B(randi(numel(c1B), numel(c1B), 1));
    i2A = c2A(randi(numel(c2A), numel(c2A), 1));
    i2B = c2B(randi(numel(c2B), numel(c2B), 1));
    w1 = svm_w(Xall, i1A, i1B, svm);
    w2 = svm_w(Xall, i2A, i2B, svm);
    if flip, w2 = -w2; end
    cosv(b) = cosine(w1, w2);
end
end

%% ------------------------------------------------------------------------
function cosv = axis_cos_relabel(Xall, c1A, c1B, c2A, c2B, nrep, svm, flip)
% Relabel null: shuffle feature labels within each condition.
cosv = zeros(nrep, 1);
i1 = [c1A; c1B]; n1A = numel(c1A);
i2 = [c2A; c2B]; n2A = numel(c2A);
for r = 1:nrep
    p1 = i1(randperm(numel(i1)));
    p2 = i2(randperm(numel(i2)));
    w1 = svm_w(Xall, p1(1:n1A), p1(n1A+1:end), svm);
    w2 = svm_w(Xall, p2(1:n2A), p2(n2A+1:end), svm);
    if flip, w2 = -w2; end
    cosv(r) = cosine(w1, w2);
end
end

%% ------------------------------------------------------------------------
function cosv = axis_cos_regroup(Xall, c1A, c1B, c2A, c2B, nrep, svm)
% Regroup null: pool by feature level across the two conditions, re-split
% into two pseudo-conditions that share one common axis by construction.
poolA = [c1A; c2A];
poolB = [c1B; c2B];
hA = floor(numel(poolA)/2);  hB = floor(numel(poolB)/2);
cosv = zeros(nrep, 1);
for r = 1:nrep
    sA = poolA(randperm(numel(poolA)));
    sB = poolB(randperm(numel(poolB)));
    w1 = svm_w(Xall, sA(1:hA),      sB(1:hB),      svm);
    w2 = svm_w(Xall, sA(hA+1:2*hA), sB(hB+1:2*hB), svm);
    cosv(r) = cosine(w1, w2);
end
end

%% ------------------------------------------------------------------------
function w = svm_w(Xall, idxA, idxB, svm)
X = [Xall(idxA, :); Xall(idxB, :)];
y = [ones(numel(idxA),1); -ones(numel(idxB),1)];
m = svm(X, y);
w = m.Beta;
end

function c = cosine(w1, w2)
c = (w1' * w2) / (norm(w1) * norm(w2) + eps);
end
