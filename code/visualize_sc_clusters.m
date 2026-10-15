function visualize_sc_clusters( ...
    features_agg, clustering, v1_classes, ...
    project_root)
% visualize_sc_clusters - Generate diagnostic figures
% for SC neuron cluster interpretation.
%
% Produces six figures:
%   1. PCA scatter (PC1 vs PC2) colored by cluster
%   2. Feature radar/bar plots per cluster
%   3. Silhouette plot
%   4. Cluster evaluation curves (k selection)
%   5. v1 vs v2 confusion matrix
%   6. Dendrogram
%
% INPUTS:
%   features_agg - Aggregated features struct from
%       extract_sc_response_features (combined across
%       sessions, quality-passing neurons only)
%       .feature_matrix  [N x F]
%       .feature_names   {1 x F}
%       .raw_features    [N x F]
%   clustering   - Output from cluster_sc_neurons
%   v1_classes   - {N x 1} cell of v1 classifications
%       ('classic', 'interneuron', 'excluded')
%   project_root - Path for saving figures
%
% Author: Claude Code
% Date: 2026-02-06

[script_dir, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'utils'));

output_dir = fullfile(project_root, 'figures');
k = clustering.k;

% Colors for clusters
if exist('richColors', 'file')
    all_colors = richColors();
else
    all_colors = lines(16);
end
cluster_colors = all_colors(1:k, :);

%% Figure 1: PCA scatter colored by cluster
fig1 = figure('Color', 'w', ...
    'Position', [100 100 900 400]);

% Panel A: colored by v2 cluster
mySubPlot([1, 2, 1]);
hold on;
for c = 1:k
    mask = clustering.assignments == c;
    scatter( ...
        clustering.pca_scores(mask, 1), ...
        clustering.pca_scores(mask, 2), ...
        30, cluster_colors(c, :), 'filled', ...
        'MarkerFaceAlpha', 0.6);
end
xlabel(sprintf('PC1 (%.1f%%)', ...
    clustering.pca_explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', ...
    clustering.pca_explained(2)));
title('v2 Clusters');
legend(clustering.labels, ...
    'Location', 'best');
hold off;

% Panel B: colored by v1 classification
mySubPlot([1, 2, 2]);
hold on;
v1_types = {'classic', 'interneuron', 'excluded'};
v1_colors = [0 0.4 0.8; 0.8 0.2 0.2; 0.5 0.5 0.5];
for i_type = 1:length(v1_types)
    mask = strcmp(v1_classes, v1_types{i_type});
    if any(mask)
        scatter( ...
            clustering.pca_scores(mask, 1), ...
            clustering.pca_scores(mask, 2), ...
            30, v1_colors(i_type, :), ...
            'filled', ...
            'MarkerFaceAlpha', 0.6);
    end
end
xlabel(sprintf('PC1 (%.1f%%)', ...
    clustering.pca_explained(1)));
ylabel(sprintf('PC2 (%.1f%%)', ...
    clustering.pca_explained(2)));
title('v1 Classes');
legend(v1_types, 'Location', 'best');
hold off;

sgtitle('PCA Projection of SC Neurons', ...
    'Interpreter', 'none');
pdfSave(fullfile(output_dir, ...
    'sc_clustering_pca.pdf'), ...
    fig1.Position(3:4) / 72, fig1);
close(fig1);

%% Figure 2: Feature bar plots per cluster
% (mean +/- SEM of raw features per cluster)
nFeatures = length(features_agg.feature_names);
n_cols = ceil(sqrt(k));
n_rows = ceil(k / n_cols);

fig2 = figure('Color', 'w', ...
    'Position', [100 100 300*n_cols 250*n_rows]);

for c = 1:k
    mySubPlot([n_rows, n_cols, c]);
    mask = clustering.assignments == c;
    n_in_cluster = sum(mask);

    % Mean and SEM of raw features
    raw_vals = features_agg.raw_features(mask, :);
    mu = mean(raw_vals, 1, 'omitnan');
    se = std(raw_vals, 0, 1, 'omitnan') / ...
        sqrt(n_in_cluster);

    bar(mu, 'FaceColor', cluster_colors(c, :), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    errorbar(1:nFeatures, mu, se, 'k.', ...
        'LineWidth', 0.5);
    hold off;

    set(gca, 'XTick', 1:nFeatures, ...
        'XTickLabel', ...
        features_agg.feature_names, ...
        'XTickLabelRotation', 60, ...
        'FontSize', 6);
    title(sprintf('%s (n=%d)', ...
        clustering.labels{c}, n_in_cluster));
    xlim([0.5, nFeatures + 0.5]);
end

sgtitle('Mean Feature Values per Cluster', ...
    'Interpreter', 'none');
pdfSave(fullfile(output_dir, ...
    'sc_clustering_features.pdf'), ...
    fig2.Position(3:4) / 72, fig2);
close(fig2);

%% Figure 3: Silhouette plot
fig3 = figure('Color', 'w', ...
    'Position', [100 100 500 400]);

[~, h] = silhouette( ...
    clustering.clust_features, ...
    clustering.assignments);
title(sprintf(['Silhouette Plot (k=%d, ' ...
    'mean=%.3f)'], k, ...
    clustering.silhouette_avg));
xlabel('Silhouette Value');
ylabel('Cluster');

pdfSave(fullfile(output_dir, ...
    'sc_clustering_silhouette.pdf'), ...
    fig3.Position(3:4) / 72, fig3);
close(fig3);

%% Figure 4: Cluster evaluation curves
fig4 = figure('Color', 'w', ...
    'Position', [100 100 900 300]);

% Panel A: Silhouette vs k
mySubPlot([1, 3, 1]);
k_range = clustering.eval_silhouette.InspectedK;
plot(k_range, ...
    clustering.eval_silhouette.CriterionValues, ...
    'ko-', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k');
hold on;
xline(clustering.eval_silhouette.OptimalK, ...
    'r--', 'LineWidth', 1);
hold off;
xlabel('k');
ylabel('Mean Silhouette');
title('Silhouette Criterion');

% Panel B: Gap vs k
mySubPlot([1, 3, 2]);
plot(k_range, ...
    clustering.eval_gap.CriterionValues, ...
    'ko-', 'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k');
hold on;
xline(clustering.eval_gap.OptimalK, ...
    'r--', 'LineWidth', 1);
hold off;
xlabel('k');
ylabel('Gap Statistic');
title('Gap Criterion');

% Panel C: BIC vs k
mySubPlot([1, 3, 3]);
bic = clustering.bic_values;
valid_k = find(isfinite(bic));
if ~isempty(valid_k)
    plot(valid_k, bic(valid_k), ...
        'ko-', 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'k');
    [~, min_idx] = min(bic(valid_k));
    hold on;
    xline(valid_k(min_idx), 'r--', ...
        'LineWidth', 1);
    hold off;
end
xlabel('k');
ylabel('BIC');
title('GMM BIC');

sgtitle(sprintf('Cluster Evaluation (selected k=%d)', ...
    k), 'Interpreter', 'none');
pdfSave(fullfile(output_dir, ...
    'sc_clustering_evaluation.pdf'), ...
    fig4.Position(3:4) / 72, fig4);
close(fig4);

%% Figure 5: v1 vs v2 confusion matrix
fig5 = figure('Color', 'w', ...
    'Position', [100 100 500 400]);

% Build crosstab
v1_unique = unique(v1_classes);
n_v1 = length(v1_unique);
conf_mat = zeros(n_v1, k);
for i = 1:n_v1
    for j = 1:k
        conf_mat(i, j) = sum( ...
            strcmp(v1_classes, v1_unique{i}) & ...
            clustering.assignments == j);
    end
end

imagesc(conf_mat);
colormap(flipud(hot));
colorbar;
set(gca, 'XTick', 1:k, ...
    'XTickLabel', clustering.labels, ...
    'XTickLabelRotation', 45, ...
    'YTick', 1:n_v1, ...
    'YTickLabel', v1_unique);
xlabel('v2 Cluster');
ylabel('v1 Class');
title('v1 vs v2 Classification');

% Add count text
for i = 1:n_v1
    for j = 1:k
        if conf_mat(i, j) > 0
            text(j, i, num2str(conf_mat(i, j)), ...
                'HorizontalAlignment', 'center', ...
                'FontWeight', 'bold', ...
                'Color', 'w');
        end
    end
end

pdfSave(fullfile(output_dir, ...
    'sc_clustering_v1_vs_v2.pdf'), ...
    fig5.Position(3:4) / 72, fig5);
close(fig5);

%% Figure 6: Dendrogram
fig6 = figure('Color', 'w', ...
    'Position', [100 100 800 400]);

% Color dendrogram by cluster assignment
dendrogram(clustering.linkage_tree, 0, ...
    'ColorThreshold', 'default');
xlabel('Neuron');
ylabel('Distance (Ward)');
title('Hierarchical Clustering Dendrogram');

pdfSave(fullfile(output_dir, ...
    'sc_clustering_dendrogram.pdf'), ...
    fig6.Position(3:4) / 72, fig6);
close(fig6);

fprintf(['visualize_sc_clusters: ' ...
    '6 figures saved to %s\n'], output_dir);

end
