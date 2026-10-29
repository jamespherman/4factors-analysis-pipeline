function res = modulation_random_time_null(spike_times, blocks, ...
    trial_events, opts)
% MODULATION_RANDOM_TIME_NULL Family-wise modulation test of one neuron
% against a random-time null (shared SC spec, section 4).
%
%   res = MODULATION_RANDOM_TIME_NULL(spike_times, blocks, ...
%       trial_events, opts)
%
%   Statistic (as in modulation_run_stat): per 25 ms test bin, the ROC
%   area of the test-bin counts across trials against the pooled
%   baseline counts; a run of n_run consecutive bins deviating in the
%   same direction scores its smallest deviation |AUC - 0.5|.
%
%   Null: n_null pseudo-blocks. In each, every trial is moved by one
%   random shift (uniform within +-opts.max_shift_s, the same shift for
%   all of the trial's events, clipped to the recording), so the
%   pseudo-trials keep the real trial count, event geometry and the
%   neuron's own burst structure, but carry no task events. Bin
%   shuffling within trials breaks burst structure and gives too narrow
%   a null for bursty neurons.
%
%   Family-wise (single-step max-T): each test (block x direction) is
%   scaled by its own null 99th percentile (floor opts.scale_floor), so
%   blocks with few trials (wide null) do not swamp the rest; the
%   largest scaled score per pseudo-block forms the null of the maximum.
%
%   INPUTS:
%       spike_times  - Spike times of the neuron (s)
%       blocks       - struct array, one per event x location:
%                      .event_idx  column of trial_events used
%                      .trials     rows of trial_events (block trials)
%                      .bin_starts bin starts relative to the event (s);
%                                  first n_baseline are baseline bins
%                      .n_baseline number of baseline bins
%                      Blocks of the same event must share bin_starts.
%       trial_events - [nTrials x nEvents] event times (s), NaN if absent
%       opts         - .bin_width, .n_run, .n_null, .max_shift_s,
%                      .rec_span [t0 t1], .event_windows [nEvents x 2]
%                      (relative windows that must fit in the
%                      recording), .scale_floor, .seed
%
%   OUTPUT res:
%       .p_block   [nBlocks x 1] family-wise p of each block
%       .dir_block [nBlocks x 1] 1 = increase, 2 = decrease
%       .p_neuron  smallest p_block
%       .run_obs   {nBlocks} [2 x nRuns] observed run scores
%       .run_p     {nBlocks} [2 x nRuns] family-wise p of each run
%       .auc       {nBlocks} [1 x nTestBins]
%
%   See also: compute_task_modulation, modulation_run_stat

n_blocks = numel(blocks);
n_trials = size(trial_events, 1);
% trials without a time for the block's event cannot be counted
for b = 1:n_blocks
    tr = blocks(b).trials(:);
    blocks(b).trials = tr(~isnan(trial_events(tr, blocks(b).event_idx)));
end
rng(opts.seed, 'twister');

% Cumulative spike count lookup, built once: number of spikes < x is
% cum(discretize(x, edges, 'IncludedEdge', 'right'))
st = sort(double(spike_times(:)));
[u, ~, j] = unique(st);
if isempty(u)
    u = 0; cum = [0; 0];
else
    cum = [0; cumsum(accumarray(j, 1))];
end
edges = [-inf; u; inf];
lut = struct('edges', edges, 'cum', cum);

% Per event: bin starts and the trials any block of that event uses
ev = unique([blocks.event_idx]);
ev_bins = cell(max(ev), 1); ev_trials = cell(max(ev), 1);
for e = ev
    b = find([blocks.event_idx] == e);
    ev_bins{e} = blocks(b(1)).bin_starts(:)';
    ev_trials{e} = unique(vertcat(blocks(b).trials));
end
% row of each block's trials within its event's trial list
block_rows = cell(n_blocks, 1);
for b = 1:n_blocks
    [~, block_rows{b}] = ismember(blocks(b).trials, ...
        ev_trials{blocks(b).event_idx});
end

%% Observed statistics
res.run_obs = cell(n_blocks, 1);
res.auc = cell(n_blocks, 1);
C_ev = event_counts(lut, ev, ev_bins, ev_trials, trial_events, ...
    zeros(n_trials, 1), opts.bin_width);
for b = 1:n_blocks
    C = C_ev{blocks(b).event_idx}(block_rows{b}, :);
    [res.run_obs{b}, ~, res.auc{b}] = modulation_run_stat(C, ...
        blocks(b).n_baseline, opts.n_run, 0);
end

%% Random-time null: one shift per trial per pseudo-block
% Allowed shift range per trial, so every event window stays inside
% the recording
lo = -opts.max_shift_s * ones(n_trials, 1);
hi = opts.max_shift_s * ones(n_trials, 1);
for e = 1:size(trial_events, 2)
    t = trial_events(:, e);
    ok = ~isnan(t);
    lo(ok) = max(lo(ok), opts.rec_span(1) - ...
        (t(ok) + opts.event_windows(e, 1)));
    hi(ok) = min(hi(ok), opts.rec_span(2) - ...
        (t(ok) + opts.event_windows(e, 2)));
end
hi = max(hi, lo);
null_score = zeros(opts.n_null, n_blocks, 2);
for k = 1:opts.n_null
    shift = lo + rand(n_trials, 1) .* (hi - lo);
    C_ev = event_counts(lut, ev, ev_bins, ev_trials, trial_events, ...
        shift, opts.bin_width);
    for b = 1:n_blocks
        C = C_ev{blocks(b).event_idx}(block_rows{b}, :);
        r = modulation_run_stat(C, blocks(b).n_baseline, opts.n_run, 0);
        null_score(k, b, :) = max(r, [], 2);
    end
end

%% Single-step max-T over blocks x directions
scale = max(prctile(null_score, 99, 1), opts.scale_floor);
null_max = max(reshape(null_score ./ scale, opts.n_null, []), [], 2);
fw_p = @(v) (1 + sum(null_max >= v)) / (1 + opts.n_null);
res.p_block = ones(n_blocks, 1);
res.dir_block = ones(n_blocks, 1);
res.run_p = cell(n_blocks, 1);
for b = 1:n_blocks
    obs = res.run_obs{b} ./ squeeze(scale(1, b, :));
    [best, res.dir_block(b)] = max(max(obs, [], 2));
    res.p_block(b) = fw_p(best);
    res.run_p{b} = arrayfun(fw_p, obs);
end
res.p_neuron = min([res.p_block; 1]);
end

function C_ev = event_counts(lut, ev, ev_bins, ev_trials, ...
    trial_events, shift, bw)
% Spike counts [nTrials x nBins] per event, each trial's event moved
% by its shift; a spike counts if strictly inside (start, start + bw),
% as in alignAndBinSpikes
C_ev = cell(max(ev), 1);
for e = ev
    tr = ev_trials{e};
    starts = trial_events(tr, e) + shift(tr) + ev_bins{e};
    C_ev{e} = n_below(lut, starts + bw, 'right') - ...
        n_below(lut, starts, 'left');
end
end

function n = n_below(lut, x, edge)
% Number of spikes < x (edge 'right') or <= x (edge 'left')
n = zeros(size(x));
n(:) = lut.cum(discretize(x(:), lut.edges, 'IncludedEdge', edge));
end
