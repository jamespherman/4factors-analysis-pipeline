function validation = validate_sc_classification( ...
    all_features, all_clustering, session_ids, ...
    v1_classes, project_root)
% validate_sc_classification - Assess quality of SC
% neuron clustering.
%
% Computes validation metrics:
%   1. Silhouette score (target > 0.4)
%   2. Bootstrap stability (target > 80%)
%   3. Cross-session stability (leave-one-out)
%   4. v1 vs v2 concordance
%
% INPUTS:
%   all_features   - Aggregated features struct with:
%       .feature_matrix  [N x F] z-scored
%       .raw_features    [N x F]
%       .feature_names   {1 x F}
%   all_clustering - Output from cluster_sc_neurons
%   session_ids    - {N x 1} cell of session IDs
%   v1_classes     - {N x 1} cell of v1 class strings
%   project_root   - Path for saving figures
%
% OUTPUT:
%   validation - Struct with metrics and summary
%
% Author: Claude Code
% Date: 2026-02-06

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

N = length(all_clustering.assignments);
k = all_clustering.k;
output_dir = fullfile(project_root, 'figures');

fprintf(['validate_sc_classification: ' ...
    'Validating %d neurons, k=%d.\n'], N, k);

%% Metric 1: Silhouette score
sil_avg = all_clustering.silhouette_avg;
sil_pass = sil_avg > 0.4;
fprintf('  Silhouette: %.3f (%s, target > 0.4)\n', ...
    sil_avg, tf_str(sil_pass));

%% Metric 2: Bootstrap stability
boot_stab = all_clustering.bootstrap_stability;
mean_boot = mean(boot_stab);
boot_pass = mean_boot > 0.80;
fprintf(['  Bootstrap stability: %.1f%% ' ...
    '(%s, target > 80%%)\n'], ...
    mean_boot * 100, tf_str(boot_pass));

%% Metric 3: Cross-session stability (leave-one-out)
unique_sessions = unique(session_ids);
n_sessions = length(unique_sessions);
loo_agreement = nan(n_sessions, 1);

if n_sessions >= 3
    fprintf(['  Cross-session leave-one-out ' ...
        '(%d sessions)...\n'], n_sessions);

    clust_features = all_clustering.clust_features;
    centroids = all_clustering.centroids;

    for i_s = 1:n_sessions
        % Identify held-out neurons
        held_out = strcmp( ...
            session_ids, unique_sessions{i_s});
        train_mask = ~held_out;

        if sum(train_mask) < 30 || ...
                sum(held_out) < 2
            continue;
        end

        % Re-cluster training set
        try
            rng(42);
            [train_assign, train_cent] = kmeans( ...
                clust_features(train_mask, :), ...
                k, ...
                'Replicates', 10, ...
                'MaxIter', 300, ...
                'Display', 'off');
        catch
            continue;
        end

        % Assign held-out neurons to nearest
        % training centroid
        held_features = ...
            clust_features(held_out, :);
        n_held = sum(held_out);
        held_assign = zeros(n_held, 1);
        for i_n = 1:n_held
            dists = sum(( ...
                train_cent - ...
                held_features(i_n, :)).^2, 2);
            [~, held_assign(i_n)] = min(dists);
        end

        % Map training clusters to original
        % assignments
        orig_train = ...
            all_clustering.assignments(train_mask);
        mapping = match_clusters_val( ...
            orig_train, train_assign, k);

        % Map held-out assignments
        mapped_held = mapping(held_assign);
        orig_held = ...
            all_clustering.assignments(held_out);

        % Compute agreement
        loo_agreement(i_s) = ...
            mean(mapped_held == orig_held);
    end

    mean_loo = mean(loo_agreement, 'omitnan');
    loo_pass = mean_loo > 0.70;
    fprintf(['  Cross-session agreement: ' ...
        '%.1f%% (%s, target > 70%%)\n'], ...
        mean_loo * 100, tf_str(loo_pass));
else
    mean_loo = NaN;
    loo_pass = false;
    fprintf(['  Cross-session: skipped ' ...
        '(only %d sessions)\n'], n_sessions);
end

%% Metric 4: v1 vs v2 concordance
% What fraction of v1 'classic' neurons are in the
% same v2 cluster?
classic_mask = strcmp(v1_classes, 'classic');
intern_mask = strcmp(v1_classes, 'interneuron');

if any(classic_mask)
    % Find dominant v2 cluster for classics
    classic_assigns = ...
        all_clustering.assignments(classic_mask);
    [mode_classic, freq_classic] = mode( ...
        classic_assigns);
    classic_purity = ...
        freq_classic / sum(classic_mask);
    fprintf(['  v1 classic purity in v2: ' ...
        '%.1f%% in cluster %d (%s)\n'], ...
        classic_purity * 100, mode_classic, ...
        all_clustering.labels{mode_classic});
else
    classic_purity = NaN;
end

if any(intern_mask)
    intern_assigns = ...
        all_clustering.assignments(intern_mask);
    [mode_intern, freq_intern] = mode( ...
        intern_assigns);
    intern_purity = ...
        freq_intern / sum(intern_mask);
    fprintf(['  v1 interneuron purity in v2: ' ...
        '%.1f%% in cluster %d (%s)\n'], ...
        intern_purity * 100, mode_intern, ...
        all_clustering.labels{mode_intern});
else
    intern_purity = NaN;
end

%% Per-cluster summary
fprintf('\n  Per-cluster summary:\n');
fprintf('  %-15s %5s %8s %8s %8s\n', ...
    'Cluster', 'N', 'Sil', 'Boot%', 'v1 mix');
for c = 1:k
    c_mask = all_clustering.assignments == c;
    c_n = sum(c_mask);
    c_sil = mean( ...
        all_clustering.silhouette_vals(c_mask));
    c_boot = mean(boot_stab(c_mask)) * 100;

    % v1 composition
    n_clas = sum(c_mask & classic_mask);
    n_int = sum(c_mask & intern_mask);
    v1_str = sprintf('%dc/%di', n_clas, n_int);

    fprintf('  %-15s %5d %8.3f %7.1f%% %8s\n', ...
        all_clustering.labels{c}, ...
        c_n, c_sil, c_boot, v1_str);
end

%% Generate validation summary figure
fig = figure('Color', 'w', ...
    'Position', [100 100 700 500]);

% Panel A: Bootstrap stability histogram
mySubPlot([2, 2, 1]);
histogram(boot_stab, 20, ...
    'FaceColor', [0.3 0.3 0.8], ...
    'EdgeColor', 'none');
xline(0.80, 'r--', 'LineWidth', 1.5);
xlabel('Bootstrap Stability');
ylabel('Count');
title(sprintf('Stability (mean=%.1f%%)', ...
    mean_boot * 100));

% Panel B: LOO agreement per session
mySubPlot([2, 2, 2]);
if ~all(isnan(loo_agreement))
    bar(loo_agreement, 'FaceColor', ...
        [0.3 0.7 0.3], 'EdgeColor', 'none');
    yline(0.70, 'r--', 'LineWidth', 1.5);
    xlabel('Session');
    ylabel('Agreement');
    title(sprintf('LOO (mean=%.1f%%)', ...
        mean_loo * 100));
    set(gca, 'XTick', 1:n_sessions, ...
        'XTickLabel', unique_sessions, ...
        'XTickLabelRotation', 45, ...
        'FontSize', 6);
else
    text(0.5, 0.5, 'Insufficient sessions', ...
        'HorizontalAlignment', 'center');
    title('LOO Agreement');
end

% Panel C: Silhouette per cluster
mySubPlot([2, 2, 3]);
cluster_sils = nan(k, 1);
for c = 1:k
    c_mask = all_clustering.assignments == c;
    cluster_sils(c) = mean( ...
        all_clustering.silhouette_vals(c_mask));
end
bar(cluster_sils, 'FaceColor', ...
    [0.8 0.5 0.2], 'EdgeColor', 'none');
yline(0.4, 'r--', 'LineWidth', 1.5);
set(gca, 'XTick', 1:k, ...
    'XTickLabel', all_clustering.labels, ...
    'XTickLabelRotation', 45);
ylabel('Mean Silhouette');
title('Per-Cluster Silhouette');

% Panel D: Summary text
mySubPlot([2, 2, 4]);
axis off;
summary_text = { ...
    sprintf('Neurons: %d', N), ...
    sprintf('Clusters: %d', k), ...
    sprintf('Silhouette: %.3f %s', ...
        sil_avg, pass_str(sil_pass)), ...
    sprintf('Bootstrap: %.1f%% %s', ...
        mean_boot * 100, pass_str(boot_pass)), ...
    sprintf('LOO: %.1f%% %s', ...
        mean_loo * 100, pass_str(loo_pass)), ...
    sprintf('Classic purity: %.1f%%', ...
        classic_purity * 100), ...
    };
text(0.1, 0.9, summary_text, ...
    'VerticalAlignment', 'top', ...
    'FontSize', 10, 'FontName', 'FixedWidth');
title('Validation Summary');

sgtitle('SC Classification Validation', ...
    'Interpreter', 'none');
pdfSave(fullfile(output_dir, ...
    'sc_classification_validation.pdf'), ...
    fig.Position(3:4) / 72, fig);
close(fig);

%% Pack output
validation = struct();
validation.silhouette_avg = sil_avg;
validation.silhouette_pass = sil_pass;
validation.bootstrap_mean = mean_boot;
validation.bootstrap_pass = boot_pass;
validation.loo_agreement = loo_agreement;
validation.loo_mean = mean_loo;
validation.loo_pass = loo_pass;
validation.classic_purity = classic_purity;
validation.intern_purity = intern_purity;

fprintf(['validate_sc_classification: ' ...
    'Validation complete.\n']);

end

%% Helper functions
function s = tf_str(tf)
if tf
    s = 'PASS';
else
    s = 'FAIL';
end
end

function s = pass_str(tf)
if tf
    s = '[PASS]';
else
    s = '[FAIL]';
end
end

function mapping = match_clusters_val( ...
    ref_assign, test_assign, k)
% Greedy cluster matching via contingency table
mapping = zeros(k, 1);
contingency = zeros(k, k);
for i = 1:k
    for j = 1:k
        contingency(i, j) = sum( ...
            ref_assign == i & test_assign == j);
    end
end

used_ref = false(k, 1);
for iter = 1:k
    [~, max_idx] = max(contingency(:));
    [ref_k, test_k] = ind2sub([k, k], max_idx);
    mapping(test_k) = ref_k;
    used_ref(ref_k) = true;
    contingency(ref_k, :) = 0;
    contingency(:, test_k) = 0;
end

unmapped = find(mapping == 0);
unused = find(~used_ref);
for i = 1:length(unmapped)
    if i <= length(unused)
        mapping(unmapped(i)) = unused(i);
    else
        mapping(unmapped(i)) = 1;
    end
end
end
