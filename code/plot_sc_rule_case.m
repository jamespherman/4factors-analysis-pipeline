function out_file = plot_sc_rule_case(session_id, cluster_id)
% PLOT_SC_RULE_CASE One-page view of an SC unit for the classic /
% non-classic rule decision (shared SC spec, rules A, B, C).
%
%   out_file = plot_sc_rule_case(session_id, cluster_id)
%
% Columns: the 4 gSac_4factors target locations, contra first. Rows:
% raster and PSTH aligned to target onset, fixation offset and saccade
% onset. Shaded windows are the screening epochs (baseline, visual,
% delay, saccade). Bottom: per location x epoch effect sizes and the
% rule A / B / C verdicts. Trials are those the screening tests use
% (gSac_4factors, inside the unit's stable period).
%
% Output: figures/sc_rule_cases/{session}_cl{cluster}_rule_case.pdf

%% Paths and data
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fullfile(findOneDrive, 'Code', ...
    '4factors-analysis-pipeline');
codes = initCodes();
S = load(fullfile(project_root, 'data', 'processed', session_id, ...
    [session_id '_session_data.mat']), 'session_data');
sd = S.session_data; clear S
ti = sd.trialInfo; ev = sd.eventTimes; c = sd.cached_friedman_results;

%% Trials, as in screen_sc_neurons
jph = find(ti.taskCode == codes.uniqueTaskCode_gSac_jph & ...
    ev.targetReillum > 0 & ev.pdsReward > 0);
f4 = find(ti.taskCode == codes.uniqueTaskCode_gSac_4factors & ...
    ~isnan(ev.targetReillum) & ev.pdsReward > 0);
all_tr = union(jph, f4);
i_cl = find(c.cluster_ids == cluster_id);
ok = ismember(all_tr, f4);
if ~isnan(c.stable_ranges(i_cl, 1))
    ok = all_tr >= c.stable_ranges(i_cl, 1) & ...
        all_tr <= c.stable_ranges(i_cl, 2);
end
F = squeeze(c.epoch_frs(i_cl, :, :))';
ok = ok(:) & all(~isnan(F), 2);
th = ti.targetTheta(all_tr) / 10;
la = c.location_angles(1:4);
if strcmp(sd.analysis.scSide, 'right')
    is_con = la > 90 & la < 270;
else
    is_con = la < 90 | la > 270;
end
% contra first, then ipsi, each by angle
[~, ord] = sortrows([~is_con(:), la(:)]);
la = la(ord); is_con = is_con(ord);
loc_grp = ord;

%% Effect sizes: epoch - baseline, per location (paired t, Bonferroni
% over 3 epochs x all location groups, as in screen_sc_neurons)
ep_names = {'visual', 'delay', 'saccade'};
alpha_corr = 0.05 / (3 * size(c.sig_by_location, 2));
D = cell(4, 1);
eff = nan(4, 3); dz = nan(4, 3); p = nan(4, 3);
for L = 1:4
    D{L} = F(ok & th(:) == la(L), 2:4) - F(ok & th(:) == la(L), 1);
    for e = 1:3
        eff(L, e) = mean(D{L}(:, e));
        dz(L, e) = eff(L, e) / std(D{L}(:, e));
        [~, p(L, e)] = ttest(D{L}(:, e));
    end
end
act = squeeze(c.sig_by_location(i_cl, loc_grp, :));
sup = squeeze(c.suppression_by_location(i_cl, loc_grp, :));

% Rule C: an ipsi activation is ignored if the best contra location's
% response in that epoch is significantly larger (Welch, one-sided)
pC = nan(4, 3);
for L = find(~is_con)
    for e = find(act(L, :))
        pC(L, e) = 1;
        for Lc = find(is_con)
            [~, pp] = ttest2(D{Lc}(:, e), D{L}(:, e), ...
                'Vartype', 'unequal', 'Tail', 'right');
            pC(L, e) = min(pC(L, e), pp);
        end
    end
end

% Rule B: preferred location / epoch (largest visual or saccade
% response) vs the opposite location, same epoch (Welch, one-sided)
R2 = eff(:, [1 3]);
[~, k] = max(R2(:));
[pL, pe] = ind2sub(size(R2), k);
pe = pe * 2 - 1;
[~, oL] = min(abs(mod(la - la(pL) - 180 + 180, 360) - 180));
[~, pB] = ttest2(D{pL}(:, pe), D{oL}(:, pe), ...
    'Vartype', 'unequal', 'Tail', 'right');

contra_act = any(act(is_con, :), 'all');
no_sup = ~any(sup, 'all');
ipsi_A = any(act(~is_con, :), 'all');
ipsi_C = any(act(~is_con, :) & ~(pC(~is_con, :) < 0.05), 'all');
verdict = @(x) ifelse(x, 'classic', 'non-classic');
vA = verdict(contra_act & ~ipsi_A & no_sup);
vC = verdict(contra_act & ~ipsi_C & no_sup);
vB = verdict(is_con(pL) & pB < 0.05 & no_sup);

%% Figure
cols = richColors();
c_con = cols(13, :);
c_ips = cols(3, :);
ep_col = [0.85 0.85 0.85; 1 0.85 0.6; 0.8 0.9 1; 0.85 1 0.85];
% alignments: event, plot window, screening epochs drawn {idx, [t0 t1]}
al = {'targetOn', [-0.2 0.4], {1, [-0.075 0.025]; 2, [0.05 0.2]}; ...
    'fixOff', [-0.4 0.2], {3, [-0.15 0.05]}; ...
    'saccadeOnset', [-0.3 0.3], {4, [-0.025 0.05]}};
bw = 0.02;
spk = sd.spikes.times(sd.spikes.clusters == cluster_id);
fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 17 11], 'Visible', 'off');
nr = 8;
ymax = zeros(1, 3);
hP = gobjects(3, 4);
for a = 1:3
    edges = al{a, 2}(1):bw:al{a, 2}(2);
    ctr = edges(1:end-1) + bw / 2;
    for L = 1:4
        tr = all_tr(ok & th(:) == la(L));
        t0 = ev.(al{a, 1})(tr);
        tr = tr(~isnan(t0)); t0 = t0(~isnan(t0));
        cc = ifelse(is_con(L), c_con, c_ips);
        % raster
        mySubPlot([nr, 4, (2 * a - 2) * 4 + L]);
        hold on
        for j = 1:numel(t0)
            s = spk(spk >= t0(j) + edges(1) & ...
                spk <= t0(j) + edges(end)) - t0(j);
            plot([s s]', repmat([j - 0.4; j + 0.4], 1, numel(s)), ...
                'Color', cc, 'LineWidth', 1);
        end
        shade(al{a, 3}, ep_col, [0 numel(t0) + 1]);
        xlim(al{a, 2}); ylim([0 numel(t0) + 1]);
        set(gca, 'XTickLabel', [], 'YTick', [], 'Layer', 'top');
        if a == 1
            title(sprintf('%.1f%s %s (n = %d)', la(L), char(176), ...
                ifelse(is_con(L), 'CONTRA', 'ipsi'), numel(t0)), ...
                'Color', cc, 'Interpreter', 'none');
        end
        if L == 1, ylabel('trials'); end
        % PSTH, mean +- SEM
        R = zeros(numel(t0), numel(ctr));
        for j = 1:numel(t0)
            R(j, :) = histcounts(spk - t0(j), edges) / bw;
        end
        m = mean(R, 1); se = std(R, 0, 1) / sqrt(max(numel(t0), 1));
        hP(a, L) = mySubPlot([nr, 4, (2 * a - 1) * 4 + L]);
        h = barStairsFill(ctr, m - se, m + se);
        set(h(2:3), 'Color', cc);
        ymax(a) = max(ymax(a), max(m + se));
        xlabel(sprintf('time from %s (s)', al{a, 1}), ...
            'Interpreter', 'none');
        if L == 1, ylabel('sp/s'); end
    end
end
for a = 1:3
    for L = 1:4
        axes(hP(a, L)); %#ok<LAXES>
        xlim(al{a, 2}); ylim([0 ymax(a) * 1.05]);
        shade(al{a, 3}, ep_col, [0 ymax(a) * 1.05]);
        set(gca, 'Layer', 'top');
    end
end

%% Table of effects
txt = {sprintf(['%s  cl %d   SC side %s   epochs: baseline ' ...
    'targetOn -75..+25 ms, visual +50..+200, delay fixOff ' ...
    '-150..+50, saccade sacOn -25..+50'], session_id, cluster_id, ...
    sd.analysis.scSide), ...
    sprintf(['Response = epoch - baseline, sp/s (dz).  * = ' ...
    'significant activation (paired t, p < %.4f, Bonferroni);  ' ...
    'v = significant suppression;  [C p] = rule C Welch p, best ' ...
    'contra > this ipsi'], alpha_corr), ''};
txt{end + 1} = sprintf('%-16s %-30s %-30s %-30s', 'location', ...
    ep_names{:});
for L = 1:4
    row = sprintf('%-16s', sprintf('%.1f %s', la(L), ...
        ifelse(is_con(L), 'contra', 'ipsi')));
    for e = 1:3
        f = ifelse(act(L, e), '*', ifelse(sup(L, e), 'v', ' '));
        cell_txt = sprintf('%+5.1f (%+.2f)%s', eff(L, e), dz(L, e), f);
        if ~isnan(pC(L, e))
            cell_txt = sprintf('%s [C p=%.2g]', cell_txt, pC(L, e));
        end
        row = [row, sprintf(' %-30s', cell_txt)]; %#ok<AGROW>
    end
    txt{end + 1} = row; %#ok<AGROW>
end
txt{end + 1} = '';
txt{end + 1} = sprintf(['Rule A (any significant ipsi activation ' ...
    '-> non-classic):  %s'], upper(vA));
txt{end + 1} = sprintf(['Rule C (ipsi activation ignored where ' ...
    'contra is significantly larger):  %s'], upper(vC));
txt{end + 1} = sprintf(['Rule B (preferred %.1f %s > opposite ' ...
    '%.1f, Welch p = %.2g):  %s'], la(pL), ep_names{pe}, la(oL), ...
    pB, upper(vB));
axes('Position', [0.05 0.01 0.9 0.2], 'Visible', 'off');
text(0, 1, txt, 'VerticalAlignment', 'top', 'FontName', 'Consolas', ...
    'FontSize', 9, 'Interpreter', 'none');

out_dir = fullfile(project_root, 'figures', 'sc_rule_cases');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
out_file = fullfile(out_dir, sprintf('%s_cl%d_rule_case.pdf', ...
    session_id, cluster_id));
exportgraphics(fig, out_file, 'ContentType', 'image', ...
    'Resolution', 200);
exportgraphics(fig, strrep(out_file, '.pdf', '.png'), 'Resolution', 110);
close(fig);
fprintf('%s\n', strjoin(txt, newline));
fprintf('Saved %s\n', out_file);
end

%% ------------------------------------------------------------------
function shade(eps, ep_col, yl)
% Draw the screening epoch windows behind the data
for k = 1:size(eps, 1)
    w = eps{k, 2};
    h = patch([w(1) w(2) w(2) w(1)], yl([1 1 2 2]), ...
        ep_col(eps{k, 1}, :), 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    uistack(h, 'bottom');
end
end

function v = ifelse(cond, a, b)
if cond, v = a; else, v = b; end
end
