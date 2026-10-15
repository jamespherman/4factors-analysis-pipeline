function plot_decoder_generalization_scatter()
% PLOT_DECODER_GENERALIZATION_SCATTER  Within-CV vs cross-condition
% accuracy, one point per session, with 95% binomial (binofit) CIs on
% both axes. 4 panels (factor x generalization dimension); colored by
% population variant (classic vs C/I pooled).
%
%   X = cross-condition generalization accuracy (both train/test
%       directions pooled into one proportion).
%   Y = within-condition 5-fold CV accuracy (both conditions pooled).
%   Both proportions share the same trial count n = nCond1 + nCond2, so
%   the CIs are directly comparable. CIs are Clopper-Pearson exact
%   (binofit). Note: CV predictions are not fully independent across
%   folds, so the Y CIs are mildly optimistic; the cross-condition X CIs
%   use disjoint train/test sets and are clean.

rng(42);
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
[~, cdefs] = define_task_conditions();
win = cdefs.window_roc_plan.epoch_windows.visual.SC;
fig_dir = fullfile(project_root, 'figures', 'decoder_axes');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

codes = initCodes();
man = readtable(fullfile(project_root, 'config', 'session_manifest.csv'));
sc = man(strcmp(man.brain_area, 'SC'), :);
nS = height(sc);

MINPC = 8; NFOLD = 5;
variants = {'classic', 'ci'};
factors  = {'salience', 'identity'};
genDims  = {'reward', 'probability'};

rc = richColors;
cVar = {rc(3,:), rc(12,:)};       % classic = blue, C/I = dark orange
vlabel = {'Classic', 'C/I pooled'};

% results{iF,iG}.(variant) = [nS x 6]: xacc xlo xhi yacc ylo yhi (NaN if skipped)
D = cell(numel(factors), numel(genDims));
for a = 1:numel(factors)
    for b = 1:numel(genDims)
        for iv = 1:numel(variants)
            D{a,b}.(variants{iv}) = nan(nS, 6);
        end
    end
end

for k = 1:nS
    sid = sc.unique_id{k};
    f = fullfile(project_root, 'data', 'processed', sid, ...
        [sid '_session_data.mat']);
    if ~exist(f, 'file'); continue; end
    S = load(f, 'session_data'); sd = S.session_data;
    nclass = {sd.analysis.screening_info.neuron_class};

    ev = sd.analysis.core_data.spikes.targetOn; tv = ev.time_vector;
    bm = (tv >= win(1)) & (tv <= win(2));
    wr = mean(ev.rates(:, :, bm), 3, 'omitnan');
    tidx = sd.analysis.core_data.trial_indices; ti = sd.trialInfo;
    loc1 = (ti.targetLocIdx(tidx) == 1);
    hr = ti.rewardDuration(tidx) > 200; lr = ti.rewardDuration(tidx) < 200;
    b1 = ti.blockNumber(tidx) == 1;     b2 = ti.blockNumber(tidx) == 2;
    bull = ti.stimType(tidx) > 2;       img = ismember(ti.stimType(tidx), [1 2]);
    lev.salience = {ti.salience(tidx) == 1, ti.salience(tidx) == 2, bull};
    lev.identity = {ti.stimType(tidx) == 1, ti.stimType(tidx) == 2, img};
    gg.reward = {hr, lr};  gg.probability = {b1, b2};

    for iv = 1:numel(variants)
        if strcmp(variants{iv}, 'classic')
            rows = strcmp(nclass, 'classic');
        else
            rows = ismember(nclass, {'classic', 'interneuron'});
        end
        X = wr(rows, :)';
        for a = 1:numel(factors)
            base = loc1 & lev.(factors{a}){3};
            yA = lev.(factors{a}){1}; yB = lev.(factors{a}){2};
            for b = 1:numel(genDims)
                g1 = gg.(genDims{b}){1}; g2 = gg.(genDims{b}){2};
                out = gen_cv_counts(X, base, yA, yB, g1, g2, MINPC, NFOLD);
                if ~isempty(out)
                    D{a,b}.(variants{iv})(k, :) = out;
                end
            end
        end
    end
    fprintf('SCAT %s done\n', sid);
end

%% Plot
fig = figure('Color','w','Position',[80 80 900 860],'Visible','off');
pj = 0;
axall = gobjects(numel(factors)*numel(genDims),1);
for a = 1:numel(factors)
    for b = 1:numel(genDims)
        pj = pj + 1;
        ax = mySubPlot([numel(factors), numel(genDims), pj]); hold on;
        axall(pj) = ax;
        plot([0 1],[0 1],'-','Color',[0.85 0.85 0.85]);      % unity
        plot([0.5 0.5],[0 1],'--','Color',[0.8 0.8 0.8]);
        plot([0 1],[0.5 0.5],'--','Color',[0.8 0.8 0.8]);
        for iv = 1:numel(variants)
            M = D{a,b}.(variants{iv});
            for k = 1:nS
                if any(isnan(M(k,:))), continue; end
                x=M(k,1); xlo=M(k,2); xhi=M(k,3);
                y=M(k,4); ylo=M(k,5); yhi=M(k,6);
                errorbar(x, y, y-ylo, yhi-y, x-xlo, xhi-x, 'o', ...
                    'Color', cVar{iv}, 'MarkerFaceColor', cVar{iv}, ...
                    'MarkerSize', 4.5, 'CapSize', 0, 'LineWidth', 0.75, ...
                    'LineStyle', 'none');
            end
        end
        axis(ax, [0.25 0.9 0.25 0.9]); axis square;
        set(ax,'TickDir','out','Color','none','XTick',0.3:0.2:0.9, ...
            'YTick',0.3:0.2:0.9);
        if b==1, ylabel(sprintf('%s\nwithin-cond CV acc', factors{a})); end
        if a==numel(factors), xlabel('cross-condition acc'); end
        if a==1, title(sprintf('generalize across %s', genDims{b}), ...
                'FontWeight','normal'); end
    end
end
% legend on panel 1
hold(axall(1),'on');
h1 = plot(axall(1), NaN, NaN, 'o', 'Color', cVar{1}, ...
    'MarkerFaceColor', cVar{1}, 'LineStyle','none');
h2 = plot(axall(1), NaN, NaN, 'o', 'Color', cVar{2}, ...
    'MarkerFaceColor', cVar{2}, 'LineStyle','none');
legend(axall(1), [h1 h2], vlabel, 'Location','southeast', ...
    'FontSize',8,'Box','off');

try; theme('light'); catch; end
pdfSave(fullfile(fig_dir,'decoder_generalization_scatter.pdf'), ...
    fig.Position(3:4)/72, fig);
close(fig);
fprintf('SCAT saved: %s\n', fullfile(fig_dir,'decoder_generalization_scatter.pdf'));
end

%% ------------------------------------------------------------------------
function out = gen_cv_counts(X, base, yA, yB, g1, g2, minpc, nfold)
% Returns [xacc xlo xhi yacc ylo yhi] or [] if too few trials.
out = [];
c1A = find(base & g1 & yA); c1B = find(base & g1 & yB);
c2A = find(base & g2 & yA); c2B = find(base & g2 & yB);
keep = @(idx) idx(~any(isnan(X(idx,:)), 2));
c1A=keep(c1A); c1B=keep(c1B); c2A=keep(c2A); c2B=keep(c2B);
if min([numel(c1A) numel(c1B) numel(c2A) numel(c2B)]) < minpc, return; end
X1 = X([c1A;c1B],:); y1 = [ones(numel(c1A),1); -ones(numel(c1B),1)];
X2 = X([c2A;c2B],:); y2 = [ones(numel(c2A),1); -ones(numel(c2B),1)];
svm = @(Xt,yt) fitcsvm(Xt,yt,'KernelFunction','linear','Standardize',true);

% cross-condition (both directions pooled)
m1 = svm(X1,y1); c12 = sum(predict(m1,X2)==y2);
m2 = svm(X2,y2); c21 = sum(predict(m2,X1)==y1);
xCorrect = c12 + c21; nTot = numel(y1) + numel(y2);

% within-condition CV (both conditions pooled)
cc1 = cv_correct(X1,y1,nfold);
cc2 = cv_correct(X2,y2,nfold);
yCorrect = cc1 + cc2;

[xacc, xci] = binofit(xCorrect, nTot);
[yacc, yci] = binofit(yCorrect, nTot);
out = [xacc xci(1) xci(2) yacc yci(1) yci(2)];
end

function nc = cv_correct(X, y, nfold)
kf = min(nfold, min(sum(y==1), sum(y==-1)));
c = cvpartition(y, 'KFold', kf);
pred = nan(size(y));
for i = 1:c.NumTestSets
    tr = training(c,i); te = test(c,i);
    m = fitcsvm(X(tr,:), y(tr), 'KernelFunction','linear','Standardize',true);
    pred(te) = predict(m, X(te,:));
end
nc = sum(pred == y);
end
