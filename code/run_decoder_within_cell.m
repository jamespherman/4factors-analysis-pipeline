function results = run_decoder_within_cell()
% RUN_DECODER_WITHIN_CELL  Within-condition decodability of salience and
% identity in each of the four reward x probability chunks at location 1.
%
%   Chunks (location 1; at loc 1 probability = block: block 1 = high
%   probability, block 2 = low probability):
%       1 HP+HR : block 1, reward > 200
%       2 HP+LR : block 1, reward < 200
%       3 LP+HR : block 2, reward > 200
%       4 LP+LR : block 2, reward < 200
%
%   For each chunk, session, and population variant (classic / C+I
%   pooled), a linear SVM decodes:
%       salience : high vs low salience  (BULLSEYE trials)
%       identity : face vs non-face      (IMAGE trials)
%   from mean target-onset firing rate (50-250 ms). Reported as 5-fold
%   cross-validated accuracy with 95% binomial (binofit) CIs.
%
%   SVM uses a UNIFORM class prior ('Prior','uniform') so the two classes
%   are weighted equally regardless of trial-count imbalance within a
%   chunk. Cells with < MINPC trials in either class are skipped (NaN).

rng(42);
[script_dir,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir,'utils'));
project_root = fileparts(script_dir);
[~,cdefs] = define_task_conditions();
win = cdefs.window_roc_plan.epoch_windows.visual.SC;

man = readtable(fullfile(project_root,'config','session_manifest.csv'));
sc = man(strcmp(man.brain_area,'SC'),:);
nS = height(sc);

MINPC = 8; NFOLD = 5;
variants = {'classic','ci'};
factors  = {'salience','identity'};
chunks   = {'HP_HR','HP_LR','LP_HR','LP_LR'};

results = struct('params', struct('window',win,'min_per_class',MINPC, ...
    'nfold',NFOLD,'prior','uniform'), 'sessions', {{}});

for k = 1:nS
    sid = sc.unique_id{k};
    f = fullfile(project_root,'data','processed',sid,[sid '_session_data.mat']);
    if ~exist(f,'file'); continue; end
    S = load(f,'session_data'); sd = S.session_data;
    nclass = {sd.analysis.screening_info.neuron_class};

    ev = sd.analysis.core_data.spikes.targetOn; tv = ev.time_vector;
    bm = (tv>=win(1)) & (tv<=win(2));
    wr = mean(ev.rates(:,:,bm),3,'omitnan');
    tidx = sd.analysis.core_data.trial_indices; ti = sd.trialInfo;
    loc1 = (ti.targetLocIdx(tidx)==1);
    hr = ti.rewardDuration(tidx)>200; lr = ti.rewardDuration(tidx)<200;
    b1 = ti.blockNumber(tidx)==1; b2 = ti.blockNumber(tidx)==2;
    bull = ti.stimType(tidx)>2; img = ismember(ti.stimType(tidx),[1 2]);
    lev.salience = {ti.salience(tidx)==1, ti.salience(tidx)==2, bull};
    lev.identity = {ti.stimType(tidx)==1, ti.stimType(tidx)==2, img};
    chunkMask = {loc1&b1&hr, loc1&b1&lr, loc1&b2&hr, loc1&b2&lr};

    sres = struct('session_id', sid);
    for iv = 1:numel(variants)
        if strcmp(variants{iv},'classic')
            rows = strcmp(nclass,'classic');
        else
            rows = ismember(nclass,{'classic','interneuron'});
        end
        X = wr(rows,:)';
        for a = 1:numel(factors)
            setmask = lev.(factors{a}){3};
            yA = lev.(factors{a}){1}; yB = lev.(factors{a}){2};
            for c = 1:numel(chunks)
                base = chunkMask{c} & setmask;
                r = cell_cv(X, base, yA, yB, MINPC, NFOLD);
                key = sprintf('%s_%s_%s', variants{iv}, factors{a}, chunks{c});
                sres.(key) = r;
            end
        end
    end
    results.sessions{end+1} = sres;
    fprintf('WCELL %s done\n', sid);
end

save(fullfile(project_root,'data','processed','decoder_within_cell.mat'), ...
    'results','-v7.3');

% Console summary: mean CV accuracy per variant/factor/chunk (sessions with data)
fprintf('\nWCELL SUMMARY (mean CV acc [n sessions ok / n with >=%d per class])\n', MINPC);
for a = 1:numel(factors)
  for iv = 1:numel(variants)
    line = sprintf('WCELL %-8s %-8s', factors{a}, variants{iv});
    for c = 1:numel(chunks)
      accs = []; nok = 0;
      for s = 1:numel(results.sessions)
        key = sprintf('%s_%s_%s', variants{iv}, factors{a}, chunks{c});
        r = results.sessions{s}.(key);
        if r.ok, accs(end+1)=r.acc; nok=nok+1; end %#ok
      end
      if isempty(accs)
        line = [line sprintf(' | %s: -- (0/%d)', chunks{c}, nS)];
      else
        line = [line sprintf(' | %s: %.2f (%d/%d)', chunks{c}, mean(accs), nok, nS)];
      end
    end
    fprintf('%s\n', line);
  end
end
end

%% ------------------------------------------------------------------------
function r = cell_cv(X, base, yA, yB, minpc, nfold)
r = struct('ok',false,'acc',NaN,'ci',[NaN NaN],'n',[NaN NaN],'reason','');
iA = find(base & yA); iB = find(base & yB);
keep = @(idx) idx(~any(isnan(X(idx,:)),2));
iA = keep(iA); iB = keep(iB);
r.n = [numel(iA) numel(iB)];
if min(r.n) < minpc, r.reason='too_few'; return; end
Xc = [X(iA,:); X(iB,:)];
y  = [ones(numel(iA),1); -ones(numel(iB),1)];
kf = min(nfold, min(numel(iA), numel(iB)));
c = cvpartition(y,'KFold',kf);
pred = nan(size(y));
for i = 1:c.NumTestSets
    tr = training(c,i); te = test(c,i);
    m = fitcsvm(Xc(tr,:), y(tr), 'KernelFunction','linear', ...
        'Standardize',true, 'Prior','uniform');
    pred(te) = predict(m, Xc(te,:));
end
nCorrect = sum(pred==y); nTot = numel(y);
[acc, ci] = binofit(nCorrect, nTot);
r.ok = true; r.acc = acc; r.ci = ci;
end
