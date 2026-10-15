function plot_decoder_generalization_axes()
% PLOT_DECODER_GENERALIZATION_AXES  Figures for the SVM generalization +
% two-null axis analysis produced by run_decoder_generalization_axes.
%
%   Figure 1: cross-condition generalization accuracy vs within-condition
%             CV vs chance, per session, for each factor x genDim, split
%             by population variant (classic / C-I pooled).
%   Figure 2: two-null axis verdict tallies (Indistinguishable /
%             Different / Inconclusive) across sessions.
%   Figure 3: two-null cosine histograms (observed / relabel / regroup)
%             for every session that yielded a verdict.

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root, 'figures', 'decoder_axes');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

L = load(fullfile(project_root, 'data', 'processed', ...
    'decoder_generalization_axes.mat'), 'results');
R = L.results;
sess = R.sessions;
nS = numel(sess);

variants = {'classic', 'ci'};
vlabel   = {'Classic', 'C/I pooled'};
factors  = {'salience', 'identity'};
genDims  = {'reward', 'probability'};

rc = richColors;
cCross = rc(3,:);    % blue
cCV    = rc(7,:);    % teal
cChance = [0.6 0.6 0.6];
verdictColors = struct('Indistinguishable', rc(8,:), ...   % green
    'Different', rc(13,:), 'Inconclusive', [0.7 0.7 0.7]); % red / gray

%% ===== FIGURE 1: generalization accuracy =====
fig1 = figure('Color','w','Position',[80 80 1000 760],'Visible','off');
pi = 0; ax1 = [];
for iF = 1:numel(factors)
    for iG = 1:numel(genDims)
        pi = pi + 1;
        ax = mySubPlot([numel(factors), numel(genDims), pi]); hold on;
        if pi == 1, ax1 = ax; end
        for iv = 1:numel(variants)
            key = sprintf('%s_%s_%s', variants{iv}, factors{iF}, genDims{iG});
            xc = iv;
            cross = []; cv = []; chance = [];
            for s = 1:nS
                if ~isfield(sess{s}, key) || ~sess{s}.(key).ok, continue; end
                g = sess{s}.(key).gen;
                cross(end+1)  = g.acc_cross_mean;        %#ok
                cv(end+1)     = mean([g.acc_cv1 g.acc_cv2]); %#ok
                chance(end+1) = g.perm_mean;             %#ok
            end
            jit = @(v,x) (x - 0.13) + 0.26*rand(size(v));
            if ~isempty(cv)
                scatter(jit(cv,xc-0.18), cv, 22, cCV, 'filled', ...
                    'MarkerFaceAlpha', 0.6);
            end
            if ~isempty(cross)
                scatter(jit(cross,xc+0.18), cross, 26, cCross, 'filled', ...
                    'MarkerFaceAlpha', 0.75);
                plot(xc+0.18, mean(cross), 'k_', 'MarkerSize', 22, 'LineWidth', 2);
            end
            if ~isempty(cv)
                plot(xc-0.18, mean(cv), 'k_', 'MarkerSize', 22, 'LineWidth', 2);
            end
        end
        yline(0.5, '--', 'Color', cChance);
        xlim([0.4 numel(variants)+0.6]); ylim([0.2 1]);
        set(ax, 'XTick', 1:numel(variants), 'XTickLabel', vlabel, 'FontSize', 8);
        if iG == 1, ylabel(sprintf('%s\naccuracy', factors{iF})); end
        if iF == 1, title(sprintf('generalize across %s', genDims{iG}), ...
                'FontWeight','normal'); end
    end
end
% legend proxy on the first panel
if ~isempty(ax1)
    hold(ax1, 'on');
    h1 = scatter(ax1, NaN, NaN, 26, cCross, 'filled');
    h2 = scatter(ax1, NaN, NaN, 22, cCV, 'filled');
    legend(ax1, [h1 h2], {'cross-condition', 'within-condition CV'}, ...
        'Location','southwest','FontSize',7,'Box','off');
end
set(findall(fig1,'Type','Axes'),'TickDir','out','Color','none');
try; theme('light'); catch; end
pdfSave(fullfile(fig_dir,'decoder_generalization_accuracy.pdf'), ...
    fig1.Position(3:4)/72, fig1);
close(fig1);

%% ===== FIGURE 2: two-null verdict tallies =====
cats = {'Indistinguishable','Different','Inconclusive'};
fig2 = figure('Color','w','Position',[80 80 1000 640],'Visible','off');
pi = 0;
for iF = 1:numel(factors)
    for iG = 1:numel(genDims)
        pi = pi + 1;
        ax = mySubPlot([numel(factors), numel(genDims), pi]); hold on;
        counts = zeros(numel(variants), numel(cats));
        for iv = 1:numel(variants)
            key = sprintf('%s_%s_%s', variants{iv}, factors{iF}, genDims{iG});
            for s = 1:nS
                if ~isfield(sess{s}, key) || ~sess{s}.(key).ok, continue; end
                ci = find(strcmp(cats, sess{s}.(key).verdict));
                if ~isempty(ci), counts(iv,ci) = counts(iv,ci) + 1; end
            end
        end
        hb = bar(counts, 'stacked', 'EdgeColor','none');
        set(hb(1),'FaceColor',verdictColors.Indistinguishable);
        set(hb(2),'FaceColor',verdictColors.Different);
        set(hb(3),'FaceColor',verdictColors.Inconclusive);
        set(ax,'XTick',1:numel(variants),'XTickLabel',vlabel,'FontSize',8);
        ylim([0 nS]);
        if iG==1, ylabel(sprintf('%s\n# sessions', factors{iF})); end
        if iF==1, title(sprintf('axis vs %s', genDims{iG}),'FontWeight','normal'); end
        if pi==1, legend(hb, cats, 'Location','northoutside', ...
                'Orientation','horizontal','FontSize',7,'Box','off'); end
    end
end
set(findall(fig2,'Type','Axes'),'TickDir','out','Color','none');
try; theme('light'); catch; end
pdfSave(fullfile(fig_dir,'decoder_axis_verdicts.pdf'), ...
    fig2.Position(3:4)/72, fig2);
close(fig2);

%% ===== FIGURE 3: per-session two-null histograms =====
edges = linspace(-1,1,25);
for iF = 1:numel(factors)
    for iG = 1:numel(genDims)
        fig3 = figure('Color','w','Position',[60 60 1100 620],'Visible','off');
        pj = 0;
        for iv = 1:numel(variants)
            for s = 1:nS
                key = sprintf('%s_%s_%s', variants{iv}, factors{iF}, genDims{iG});
                pj = pj + 1;
                ax = mySubPlot([numel(variants), nS, pj], ...
                    'HeightSpacing', 0.11, 'WidthSpacing', 0.02); hold on;
                if ~isfield(sess{s}, key) || ~sess{s}.(key).ok
                    text(0.5,0.5,'n/a','Units','normalized', ...
                        'HorizontalAlignment','center','Color',[.6 .6 .6]);
                    set(ax,'XTick',[],'YTick',[]); continue;
                end
                a = sess{s}.(key).axis;
                hn = @(x,c,al) histogram(x,edges,'Normalization','probability', ...
                    'FaceColor',c,'EdgeColor','none','FaceAlpha',al);
                hn(a.relabelCos,[0.7 0.7 0.7],0.7);
                hn(a.regroupCos,rc(12,:),0.5);
                hn(a.obsCos,rc(3,:),0.7);
                xlim([-1 1]);
                title(sprintf_short(sess{s}.session_id, sess{s}.(key).verdict), ...
                    'FontSize',6,'FontWeight','normal','Interpreter','none');
                if s==1, ylabel(vlabel{iv},'FontSize',7); end
                set(ax,'FontSize',6);
            end
        end
        sgtitle(sprintf('%s axis across %s (blue=obs, orange=regroup, gray=relabel)', ...
            factors{iF}, genDims{iG}), 'FontSize', 9);
        set(findall(fig3,'Type','Axes'),'TickDir','out','Color','none');
        try; theme('light'); catch; end
        pdfSave(fullfile(fig_dir, sprintf('decoder_axis_hist_%s_%s.pdf', ...
            factors{iF}, genDims{iG})), fig3.Position(3:4)/72, fig3);
        close(fig3);
    end
end

fprintf('PLOTS saved to %s\n', fig_dir);
end

function s = sprintf_short(sid, verdict)
% e.g. 'Feynman_08_12_2025_SC' -> 'F 08-12 / Indistinguishable'
parts = strsplit(sid, '_');
s = sprintf('%s %s-%s\n%s', parts{1}(1), parts{2}, parts{3}, verdict);
end
