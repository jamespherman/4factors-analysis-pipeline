function [run_obs, null_max, auc] = modulation_run_stat(counts, ...
    n_baseline, n_run, n_perm)
% MODULATION_RUN_STAT Run-of-bins modulation statistic and its
% within-trial permutation null.
%
%   [run_obs, null_max, auc] = MODULATION_RUN_STAT(counts, n_baseline,
%       n_run, n_perm)
%
%   The statistic is the one used by compute_task_modulation: for each
%   test bin, the ROC area of the test-bin spike counts (across trials)
%   against the baseline counts pooled over trials and baseline bins;
%   AUC > 0.5 means more spikes in the test bin. A run is n_run
%   consecutive test bins that all deviate in the same direction, and
%   its strength is the smallest deviation |AUC - 0.5| in the run.
%
%   The null comes from permuting the bins within each trial. When the
%   neuron is not modulated, the baseline and test bins of a trial are
%   exchangeable whatever the trial's overall rate, so the permutations
%   reproduce the statistic's chance distribution for this neuron's
%   own spike counts and trial-to-trial rate changes. This matters at
%   low rates, where most bins hold 0 spikes and the ROC area can
%   differ from 0.5 by chance in every bin at once.
%
%   INPUTS:
%       counts     - [nTrials x nBins] spike counts; the first
%                    n_baseline columns are the baseline window, the
%                    rest are the test window
%       n_baseline - Number of baseline bins
%       n_run      - Consecutive bins required (e.g. 3)
%       n_perm     - Number of permutations (0 = no null)
%
%   OUTPUTS:
%       run_obs  - [2 x nRuns] run strength for each run start; row 1
%                  increases (min of AUC - 0.5), row 2 decreases
%                  (min of 0.5 - AUC). Negative = no consistent run.
%       null_max - [n_perm x 1] largest run strength (either
%                  direction, any run start) in each permutation
%       auc      - [1 x nTestBins] ROC area of each test bin
%
%   See also: compute_task_modulation

counts = double(counts);
[n_trials, n_bins] = size(counts);
n_test = n_bins - n_baseline;

auc = bin_auc(counts, n_baseline);
run_obs = run_strength(auc, n_run);

null_max = zeros(n_perm, 1);
if n_perm < 1 || n_test < n_run
    return;
end

% Permute in chunks to bound memory (trials x bins x perms)
chunk = max(1, floor(4e6 / (n_trials * n_bins)));
for p0 = 1:chunk:n_perm
    np = min(chunk, n_perm - p0 + 1);
    [~, order] = sort(rand(n_trials, n_bins, np), 2);
    rows = repmat((1:n_trials)', 1, n_bins, np);
    perm_counts = counts(rows + (order - 1) * n_trials);
    auc_p = bin_auc(perm_counts, n_baseline);
    for k = 1:np
        r = run_strength(auc_p(:, :, k), n_run);
        null_max(p0 + k - 1) = max(r(:));
    end
end

end

function auc = bin_auc(counts, n_baseline)
% ROC area of each test bin vs the pooled baseline, from the count
% histograms (counts are small integers). counts may be a stack of
% [trials x bins x perms]; auc is then [1 x nTest x nPerm].
n_trials = size(counts, 1);
n_bins = size(counts, 2);
n_stack = size(counts, 3);
n_test = n_bins - n_baseline;
base = counts(:, 1:n_baseline, :);
test = counts(:, (n_baseline + 1):end, :);
max_count = max(counts(:));
auc = zeros(1, n_test, n_stack);
% fraction of baseline counts strictly below the current value
below = zeros(1, 1, n_stack);
for v = 0:max_count
    p_base = sum(sum(base == v, 1), 2) / (n_trials * n_baseline);
    p_test = sum(test == v, 1) / n_trials;
    auc = auc + p_test .* (below + 0.5 * p_base);
    below = below + p_base;
end
end

function r = run_strength(auc, n_run)
% Smallest same-direction deviation within each run of n_run bins
dev = auc(:)' - 0.5;
n_runs = numel(dev) - n_run + 1;
if n_runs < 1
    r = -inf(2, 1);
    return;
end
idx = (1:n_runs)' + (0:(n_run - 1));
r = [min(dev(idx), [], 2)'; min(-dev(idx), [], 2)'];
end
