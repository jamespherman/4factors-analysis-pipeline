function plot_decoder_within_cell()
% PLOT_DECODER_WITHIN_CELL  Per-session within-cell decodability across the
% four reward x probability chunks at location 1. 4 panels
% (feature x population variant); x = chunk, y = 5-fold CV accuracy with
% 95% binomial (binofit) CIs; one point per session, colored by monkey.

[script_dir,~,~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir,'utils'));
project_root = fileparts(script_dir);
fig_dir = fullfile(project_root,'figures','decoder_axes');
if ~exist(fig_dir,'dir'); mkdir(fig_dir); end

L = load(fullfile(project_root,'data','processed','decoder_within_cell.mat'),'results');
R = L.results; sess = R.sessions; nS = numel(sess);

variants = {'classic','ci'}; vlabel = {'Classic','C/I pooled'};
factors  = {'salience','identity'};
chunks   = {'HP_HR','HP_LR','LP_HR','LP_LR'};
clabel   = {'HP+HR','HP+LR','LP+HR','LP+LR'};

rc = richColors;
monkeyColor = containers.Map({'Feynman','Newton'}, {rc(3,:), rc(12,:)});

fig = figure('Color','w','Position',[80 80 980 780],'Visible','off');
pj = 0; axall = gobjects(4,1);
for a = 1:numel(factors)
    for iv = 1:numel(variants)
        pj = pj + 1;
        ax = mySubPlot([numel(factors), numel(variants), pj]); hold on;
        axall(pj) = ax;
        yline(0.5,'--','Color',[0.8 0.8 0.8]);
        for c = 1:numel(chunks)
            accs = [];
            for s = 1:nS
                key = sprintf('%s_%s_%s', variants{iv}, factors{a}, chunks{c});
                r = sess{s}.(key);
                if ~r.ok, continue; end
                mk = strtok(sess{s}.session_id, '_');
                col = monkeyColor(mk);
                xj = c + 0.12*(rand-0.5)*2;
                errorbar(xj, r.acc, r.acc-r.ci(1), r.ci(2)-r.acc, 'o', ...
                    'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 4, ...
                    'CapSize', 0, 'LineWidth', 0.6, 'LineStyle','none');
                accs(end+1) = r.acc; %#ok
            end
            % mean line per chunk
            if ~isempty(accs)
                plot(c+[-0.28 0.28], [mean(accs) mean(accs)], 'k-', ...
                    'LineWidth', 1.5);
                text(c, 0.955, sprintf('%d/%d', numel(accs), nS), ...
                    'HorizontalAlignment','center','FontSize',6,'Color',[.4 .4 .4]);
            end
        end
        xlim([0.5 4.5]); ylim([0.25 1]);
        set(ax,'XTick',1:4,'XTickLabel',clabel,'FontSize',8,'TickDir','out','Color','none');
        if iv==1, ylabel(sprintf('%s\nCV accuracy', factors{a})); end
        title(vlabel{iv},'FontWeight','normal');
    end
end
% monkey legend
hold(axall(1),'on');
hF = plot(axall(1),NaN,NaN,'o','Color',monkeyColor('Feynman'), ...
    'MarkerFaceColor',monkeyColor('Feynman'),'LineStyle','none');
hN = plot(axall(1),NaN,NaN,'o','Color',monkeyColor('Newton'), ...
    'MarkerFaceColor',monkeyColor('Newton'),'LineStyle','none');
legend(axall(1),[hF hN],{'Feynman','Newton'},'Location','southwest', ...
    'FontSize',7,'Box','off');

try; theme('light'); catch; end
pdfSave(fullfile(fig_dir,'decoder_within_cell_accuracy.pdf'), ...
    fig.Position(3:4)/72, fig);
close(fig);
fprintf('WCELLPLOT saved\n');
end
