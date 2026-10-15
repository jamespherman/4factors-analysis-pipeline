function generate_factor_psths(varargin)
% GENERATE_FACTOR_PSTHS Per-neuron PSTHs for window-ROC significant factors.
%
%   generate_factor_psths()
%   generate_factor_psths('SessionID', 'Feynman_08_12_2025_SC')
%   generate_factor_psths('ForceRerun', true)
%
%   Produces one figure per SC neuron that is classified 'classic' or
%   'interneuron' AND is significantly modulated by at least one factor
%   in the window-ROC visual epoch. Significance is defined as the
%   bootstrap CI excluding 0.5, matching window_roc_scatter_SC.pdf.
%
%   LAYOUT
%   Columns are the four target locations. Each significant factor
%   contributes two rows: a raster (trials grouped by factor level and
%   colored by it) and a PSTH (one trace per level, gray fill between).
%   All activity is aligned on target onset; the shaded band marks the
%   window in which the window ROC was computed.
%
%   Probability is only contrasted at locations 1 and 3. The design makes
%   high probability occur at location 1 (block 1) and location 3
%   (block 2) only, so every trial at locations 2 and 4 is low
%   probability and no contrast exists there. Those panels are left
%   blank and annotated.
%
%   The column matching the neuron's preferred location is boxed: that is
%   the only location whose trials entered the ROC that selected this
%   neuron.
%
%   OPTIONAL INPUTS (Name-Value pairs):
%       'SessionID'  - Process only this session (default: all SC)
%       'ForceRerun' - Regenerate existing figures (default: false)
%
%   OUTPUT:
%       figures/factor_psths/{session_id}_{class}_{cluster}_factors.pdf

%% Parse inputs
p = inputParser;
addParameter(p, 'SessionID', '', @ischar);
addParameter(p, 'ForceRerun', false, @islogical);
parse(p, varargin{:});
session_filter = p.Results.SessionID;
force_rerun = p.Results.ForceRerun;

%% Setup paths
[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));
project_root = fileparts(script_dir);

output_dir = fullfile(project_root, 'figures', 'factor_psths');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Analysis configuration
[~, condition_defs] = define_task_conditions();

% The ROC window we are illustrating (visual epoch, SC)
roc_window = condition_defs.window_roc_plan.epoch_windows.visual.SC;

% Display window for the target-aligned PSTHs
t_start = -0.25;
t_end   = 0.50;
bin_width = 0.025;

% Factors, in display order. cond1 is "low", cond2 is "high", matching
% define_task_conditions: AUC > 0.5 means higher firing for "high".
factor_order = {'reward', 'salience', 'probability', 'identity'};
factor_levels = struct(...
    'reward',      {{'low reward',  'high reward'}}, ...
    'salience',    {{'low salience', 'high salience'}}, ...
    'probability', {{'low prob.',   'high prob.'}}, ...
    'identity',    {{'non-face',    'face'}});

% Probability is only manipulated at these locations
prob_locations = [1, 3];

%% Colors
colors = richColors;
col_low  = colors(3, :);    % Royal Blue
col_high = colors(13, :);   % Bright Red

%% Load manifest
manifest = readtable(fullfile(project_root, 'config', ...
    'session_manifest.csv'));
sc_sessions = manifest(strcmp(manifest.brain_area, 'SC'), :);

if ~isempty(session_filter)
    sc_sessions = sc_sessions(strcmp(sc_sessions.unique_id, ...
        session_filter), :);
    if isempty(sc_sessions)
        error('generate_factor_psths:sessionNotFound', ...
            'Session %s not found or is not an SC session.', ...
            session_filter);
    end
end

fprintf('Processing %d SC session(s)...\n', height(sc_sessions));

codes = initCodes();
total_neurons = 0;
total_figures = 0;

%% Process each session
for i_session = 1:height(sc_sessions)
    session_id = sc_sessions.unique_id{i_session};
    fprintf('\n=== %s ===\n', session_id);

    data_path = fullfile(project_root, 'data', 'processed', ...
        session_id, [session_id, '_session_data.mat']);
    if ~exist(data_path, 'file')
        fprintf('  Session data not found, skipping.\n');
        continue;
    end

    load(data_path, 'session_data');

    if ~isfield(session_data, 'analysis') || ...
       ~isfield(session_data.analysis, 'window_roc') || ...
       ~isfield(session_data.analysis, 'screening_info')
        fprintf('  Missing window_roc or screening_info, skipping.\n');
        continue;
    end

    screening_info = session_data.analysis.screening_info;
    window_roc = session_data.analysis.window_roc;

    if ~isfield(window_roc, 'visual')
        fprintf('  No visual epoch in window_roc, skipping.\n');
        continue;
    end
    wroc = window_roc.visual;

    % Row i of window_roc corresponds to cluster_info row i, which is
    % screening_info(i). This alignment depends on prepare_neuronal_data
    % indexing by cluster_info.cluster_id rather than by mask position.
    cluster_ids = session_data.spikes.cluster_info.cluster_id;
    n_neurons = numel(cluster_ids);
    if numel(screening_info) ~= n_neurons
        fprintf(['  screening_info (%d) does not match cluster count ' ...
            '(%d), skipping.\n'], numel(screening_info), n_neurons);
        continue;
    end
    if numel(wroc.(factor_order{1}).auc) ~= n_neurons
        fprintf(['  window_roc rows (%d) do not match cluster count ' ...
            '(%d), skipping.\n'], ...
            numel(wroc.(factor_order{1}).auc), n_neurons);
        continue;
    end

    %% Trial selection
    trialInfo = session_data.trialInfo;
    eventTimes = session_data.eventTimes;

    master_mask = (trialInfo.taskCode == ...
        codes.uniqueTaskCode_gSac_4factors) & ...
        ~cellfun(@isempty, eventTimes.rewardCell);
    trial_indices = find(master_mask);

    if isempty(trial_indices)
        fprintf('  No valid 4factors trials, skipping.\n');
        continue;
    end

    target_on = eventTimes.targetOn(trial_indices);
    target_loc = trialInfo.targetLocIdx(trial_indices);
    block_num = trialInfo.blockNumber(trial_indices);

    % Factor level masks, matching define_task_conditions exactly
    lvl = struct();
    lvl.reward = {trialInfo.rewardDuration(trial_indices) < 200, ...
                  trialInfo.rewardDuration(trial_indices) > 200};

    is_bullseye = trialInfo.stimType(trial_indices) > 2;
    lvl.salience = {...
        (trialInfo.salience(trial_indices) == 2) & is_bullseye, ...
        (trialInfo.salience(trial_indices) == 1) & is_bullseye};

    is_image = ismember(trialInfo.stimType(trial_indices), [1, 2]);
    lvl.identity = {...
        (trialInfo.stimType(trial_indices) == 2) & is_image, ...
        (trialInfo.stimType(trial_indices) == 1) & is_image};

    is_high_prob = (block_num == 1 & target_loc == 1) | ...
                   (block_num == 2 & target_loc == 3);
    is_low_prob = (block_num == 1 & ismember(target_loc, [2, 3, 4])) | ...
                  (block_num == 2 & ismember(target_loc, [1, 2, 4]));
    lvl.probability = {is_low_prob, is_high_prob};

    % Per-neuron stable trial masks, matching analyze_window_roc
    if isfield(session_data, 'stability_info') && ...
            ~isempty(session_data.stability_info)
        neuron_trial_masks = get_neuron_present_trials(...
            session_data.stability_info, trial_indices, n_neurons);
    else
        neuron_trial_masks = true(n_neurons, numel(trial_indices));
    end

    %% Identify neurons to plot
    neuron_class = {screening_info.neuron_class};
    is_candidate = ismember(neuron_class, {'classic', 'interneuron'}) & ...
        [screening_info.included];

    % Significance: bootstrap CI excludes 0.5 (matches the scatter figure)
    sig_by_factor = false(n_neurons, numel(factor_order));
    for i_f = 1:numel(factor_order)
        f = factor_order{i_f};
        if ~isfield(wroc, f), continue; end
        ci = wroc.(f).ci;
        s = (ci(:, 1) > 0.5) | (ci(:, 2) < 0.5);
        s(isnan(ci(:, 1)) | isnan(ci(:, 2))) = false;
        sig_by_factor(:, i_f) = s;
    end

    to_plot = find(is_candidate(:) & any(sig_by_factor, 2));
    fprintf('  %d/%d neurons: C/I and significant for >=1 factor.\n', ...
        numel(to_plot), n_neurons);
    total_neurons = total_neurons + numel(to_plot);

    %% Build a figure per neuron
    for i_sel = 1:numel(to_plot)
        i_neuron = to_plot(i_sel);
        cluster_id = cluster_ids(i_neuron);
        nclass = neuron_class{i_neuron};

        fig_name = sprintf('%s_%s_%03d_factors.pdf', session_id, ...
            nclass, cluster_id);
        fig_path = fullfile(output_dir, fig_name);
        if exist(fig_path, 'file') && ~force_rerun
            continue;
        end

        spike_times = session_data.spikes.times(...
            session_data.spikes.clusters == cluster_id);
        if isempty(spike_times)
            fprintf('    cluster %d: no spikes, skipping.\n', cluster_id);
            continue;
        end

        stable_mask = neuron_trial_masks(i_neuron, :)';
        sig_factors = factor_order(sig_by_factor(i_neuron, :));

        ok = plot_neuron_figure(fig_path, spike_times, target_on, ...
            target_loc, stable_mask, lvl, sig_factors, factor_levels, ...
            wroc, i_neuron, session_id, cluster_id, nclass, ...
            col_low, col_high, t_start, t_end, bin_width, ...
            roc_window, prob_locations);

        if ok
            total_figures = total_figures + 1;
            fprintf('    Saved %s (%s)\n', fig_name, ...
                strjoin(sig_factors, ', '));
        end
    end
end

fprintf('\n=== Complete ===\n');
fprintf('Neurons meeting criteria: %d\n', total_neurons);
fprintf('Figures generated: %d\n', total_figures);
fprintf('Output: %s\n', output_dir);

end


%% ========================================================================
%  PER-NEURON FIGURE
%  ========================================================================
function ok = plot_neuron_figure(fig_path, spike_times, target_on, ...
    target_loc, stable_mask, lvl, sig_factors, factor_levels, ...
    wroc, i_neuron, session_id, cluster_id, nclass, ...
    col_low, col_high, t_start, t_end, bin_width, roc_window, ...
    prob_locations)

ok = false;
n_fac = numel(sig_factors);
n_rows = 2 * n_fac;
n_cols = 4;

% Reserve fixed pixel bands for the header and the x-axis label rather
% than relying on mySubPlot's proportional defaults, which leave too
% little room on short (single-factor) figures.
panel_h = 150;
top_reserve = 46;
bot_reserve = 52;
fig_h = panel_h * n_rows + top_reserve + bot_reserve;
fig_w = 240 * n_cols;
sub_opts = {'Height', (panel_h * n_rows) / fig_h, ...
            'BottomMargin', bot_reserve / fig_h};

fig = figure('Color', 'w', 'Visible', 'off', ...
    'Position', [100, 100, fig_w, fig_h]);

% First pass: gather firing rates so PSTH y-limits can be shared within
% each factor row.
psth_store = cell(n_fac, n_cols, 2);
rast_store = cell(n_fac, n_cols, 2);
time_vec = [];

for i_f = 1:n_fac
    f = sig_factors{i_f};
    for loc = 1:n_cols
        if strcmp(f, 'probability') && ~ismember(loc, prob_locations)
            continue;
        end
        loc_mask = (target_loc == loc) & stable_mask;
        for i_lev = 1:2
            trial_mask = loc_mask & lvl.(f){i_lev};
            ev = target_on(trial_mask);
            ev = ev(~isnan(ev));
            rast_store{i_f, loc, i_lev} = align_spikes(spike_times, ...
                ev, t_start, t_end);
            [tv, fr] = binned_rate(spike_times, ev, t_start, ...
                t_end, bin_width);
            % A level with no trials returns an empty time vector; do
            % not let it clobber a good one from another panel.
            if ~isempty(tv)
                time_vec = tv;
            end
            psth_store{i_f, loc, i_lev} = fr;
        end
    end
end

if isempty(time_vec)
    close(fig);
    return;
end

%% Draw
for i_f = 1:n_fac
    f = sig_factors{i_f};
    levels = factor_levels.(f);
    pref_loc = get_preferred_location(wroc, f, i_neuron);

    % Shared limits for this factor's rows
    y_max = 0;
    for loc = 1:n_cols
        for i_lev = 1:2
            fr = psth_store{i_f, loc, i_lev};
            if ~isempty(fr), y_max = max(y_max, max(fr)); end
        end
    end
    if y_max <= 0, y_max = 1; end
    y_max = ceil(y_max / 5) * 5;

    n_trials_max = 1;
    for loc = 1:n_cols
        n_here = 0;
        for i_lev = 1:2
            n_here = n_here + numel(rast_store{i_f, loc, i_lev});
        end
        n_trials_max = max(n_trials_max, n_here);
    end

    for loc = 1:n_cols
        row_rast = (i_f - 1) * 2 + 1;
        row_psth = (i_f - 1) * 2 + 2;

        % --- Raster ---
        ax_r = mySubPlot([n_rows, n_cols, (row_rast - 1) * n_cols + loc], ...
            sub_opts{:});
        hold on;

        if strcmp(f, 'probability') && ~ismember(loc, prob_locations)
            draw_na_panel(ax_r, t_start, t_end);
        else
            y0 = 0;
            for i_lev = 1:2
                trials = rast_store{i_f, loc, i_lev};
                c = col_low;
                if i_lev == 2, c = col_high; end
                draw_raster(trials, y0, c);
                y0 = y0 + numel(trials);
            end
            shade_window(roc_window, [0, max(n_trials_max, 1)]);
            plot([0 0], [0, max(n_trials_max, 1)], 'k--', ...
                'LineWidth', 0.5);
            ylim([0, max(n_trials_max, 1)]);
            xlim([t_start, t_end]);
            set(ax_r, 'XTickLabel', []);
            if loc == 1
                ylabel('trials');
            else
                % Panels abut, so repeated tick labels collide with the
                % neighbouring axes. Scale is shared across the row.
                set(ax_r, 'YTickLabel', []);
            end
        end

        % Column headers on the very top row
        if i_f == 1
            title(sprintf('Location %d', loc), 'FontWeight', 'normal');
        end

        % Box the location that actually entered the ROC
        if loc == pref_loc
            set(ax_r, 'LineWidth', 2, 'XColor', [0.2 0.2 0.2], ...
                'YColor', [0.2 0.2 0.2], 'Box', 'on');
        end

        % --- PSTH ---
        ax_p = mySubPlot([n_rows, n_cols, (row_psth - 1) * n_cols + loc], ...
            sub_opts{:});
        hold on;

        if strcmp(f, 'probability') && ~ismember(loc, prob_locations)
            draw_na_panel(ax_p, t_start, t_end);
        else
            fr_low = psth_store{i_f, loc, 1};
            fr_high = psth_store{i_f, loc, 2};
            shade_window(roc_window, [0, y_max]);

            has_low = ~isempty(fr_low) && ~all(isnan(fr_low));
            has_high = ~isempty(fr_high) && ~all(isnan(fr_high));
            if has_low && has_high
                % Both levels present: fill the gap between them so the
                % size of the effect reads at a glance.
                h = barStairsFill(time_vec, fr_low, fr_high);
                set(h(1), 'FaceColor', 0.85 * [1 1 1]);
                set(h(2), 'Color', col_low, 'LineWidth', 1.25);
                set(h(3), 'Color', col_high, 'LineWidth', 1.25);
            elseif has_low
                set(barStairs(time_vec, fr_low, false), ...
                    'Color', col_low, 'LineWidth', 1.25);
            elseif has_high
                set(barStairs(time_vec, fr_high, false), ...
                    'Color', col_high, 'LineWidth', 1.25);
            end

            plot([0 0], [0, y_max], 'k--', 'LineWidth', 0.5);
            ylim([0, y_max]);
            xlim([t_start, t_end]);
            if loc == 1
                ylabel(sprintf('%s\nsp/s', f), 'Interpreter', 'none');
                % Level names differ by factor, so give each its own key.
                h1 = plot(NaN, NaN, 'Color', col_low, 'LineWidth', 1.25);
                h2 = plot(NaN, NaN, 'Color', col_high, 'LineWidth', 1.25);
                legend([h1, h2], levels, 'Location', 'northwest', ...
                    'FontSize', 6, 'Box', 'off');
            else
                set(ax_p, 'YTickLabel', []);
            end
            if loc == pref_loc
                auc = wroc.(f).auc(i_neuron);
                ci = wroc.(f).ci(i_neuron, :);
                title(sprintf('AUC=%.2f [%.2f %.2f]', auc, ci(1), ...
                    ci(2)), 'FontWeight', 'normal', 'FontSize', 8);
            end
            if i_f == n_fac
                xlabel('time from target onset (s)');
            end
        end

        if loc == pref_loc
            set(ax_p, 'LineWidth', 2, 'XColor', [0.2 0.2 0.2], ...
                'YColor', [0.2 0.2 0.2], 'Box', 'on');
        end
    end
end

%% Header
annotation(fig, 'textbox', ...
    [0.01, 1 - 30 / fig_h, 0.98, 24 / fig_h], 'String', ...
    sprintf(['%s  |  cluster %d  |  %s  |  target-aligned  |  ' ...
    'shaded = window ROC epoch [%.0f-%.0f ms]  |  ' ...
    'boxed column = preferred location (the only one in the ROC)'], ...
    session_id, cluster_id, nclass, roc_window(1) * 1000, ...
    roc_window(2) * 1000), ...
    'Interpreter', 'none', 'EdgeColor', 'none', 'FontSize', 9, ...
    'HorizontalAlignment', 'left');

set(findall(fig, 'Type', 'Axes'), 'TickDir', 'out', 'Color', 'none');

try
    theme('light');
catch
end

pdfSave(fig_path, fig.Position(3:4) / 72, fig);
close(fig);
ok = true;

end


%% ========================================================================
%  HELPERS
%  ========================================================================
function pref = get_preferred_location(wroc, factor_name, i_neuron)
% Location whose trials entered the window ROC for this neuron.
pref = NaN;
if isfield(wroc, factor_name) && ...
        isfield(wroc.(factor_name), 'preferred_location')
    pref = wroc.(factor_name).preferred_location(i_neuron);
end
end


function trials = align_spikes(spike_times, event_times, t_start, t_end)
% Per-trial spike times relative to each event, for raster plotting.
trials = cell(numel(event_times), 1);
for i = 1:numel(event_times)
    t0 = event_times(i);
    s = spike_times(spike_times >= t0 + t_start & ...
        spike_times <= t0 + t_end) - t0;
    trials{i} = s(:)';
end
end


function [time_vec, fr] = binned_rate(spike_times, event_times, ...
    t_start, t_end, bin_width)
% Mean firing rate in non-overlapping bins, aligned to event_times.
time_vec = [];
fr = [];
if isempty(event_times)
    return;
end
[counts, centers, ~] = alignAndBinSpikes(spike_times, event_times, ...
    t_start, t_end, bin_width);
fr = mean(counts, 1, 'omitnan') / bin_width;
fr = fr(:)';
time_vec = centers(:)';
end


function draw_raster(trials, y_offset, color)
% Vertical tick per spike; one row per trial.
if isempty(trials), return; end
xs = [];
ys = [];
for r = 1:numel(trials)
    s = trials{r};
    if isempty(s), continue; end
    n = numel(s);
    xs = [xs, [s; s; nan(1, n)]];
    ys = [ys, [repmat(y_offset + r - 0.45, 1, n); ...
               repmat(y_offset + r + 0.45, 1, n); ...
               nan(1, n)]];
end
if isempty(xs), return; end
plot(xs(:), ys(:), 'Color', color, 'LineWidth', 0.5);
end


function shade_window(win, ylims)
% Shade the ROC analysis window behind the data.
h = fill([win(1), win(2), win(2), win(1)], ...
    [ylims(1), ylims(1), ylims(2), ylims(2)], ...
    [0.95, 0.92, 0.70], 'EdgeColor', 'none');
uistack(h, 'bottom');
end


function draw_na_panel(ax, t_start, t_end)
% Probability is not contrasted at this location; say so rather than
% implying an absent effect.
xlim(ax, [t_start, t_end]);
ylim(ax, [0, 1]);
text(0.5, 0.5, {'high prob.', 'never occurs', 'here'}, ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontSize', 7, ...
    'Color', [0.45 0.45 0.45]);
set(ax, 'XTick', [], 'YTick', [], 'Box', 'off', ...
    'XColor', 'none', 'YColor', 'none');
end
