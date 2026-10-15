function features = extract_sc_response_features( ...
    session_data, scSide)
% extract_sc_response_features - Extract per-location
% response features for SC neuron clustering.
%
% Uses cached_friedman_results to build a feature vector
% for each neuron. Locations are rank-ordered by visual
% response strength (best to worst) so that features are
% RF-position invariant across neurons and sessions.
%
% Only the 4 gSac_4factors locations are used (the 5th
% gSac_jph group, if present, is excluded).
%
% INPUTS:
%   session_data - Session data struct with fields:
%       .cached_friedman_results
%       .spikes.cluster_info.cluster_id
%   scSide - String ('right' or 'left'), passed through
%            but not used for feature computation.
%
% OUTPUT:
%   features - Struct with fields:
%       .cluster_ids     [N x 1] cluster IDs
%       .feature_matrix  [N x F] z-scored features
%       .raw_features    [N x F] before normalization
%       .feature_names   {1 x F} cell of feature names
%       .passes_quality  [N x 1] logical
%       .scSide          string, passed through
%
% Author: Claude Code
% Date: 2026-02-06

%% Validate inputs
if ~isfield(session_data, ...
        'cached_friedman_results') ...
        || isempty( ...
        session_data.cached_friedman_results)
    error(['extract_sc_response_features: ' ...
        'No cached_friedman_results found.']);
end

cached = session_data.cached_friedman_results;
cluster_ids = cached.cluster_ids;
nClusters = length(cluster_ids);

%% Extract data from cached results
location_angles = cached.location_angles;
n_groups = length(location_angles);
group_mean_frs = cached.group_mean_frs;
sig_by_loc = cached.sig_by_location;
supp_by_loc = cached.suppression_by_location;

% Use only first 4 groups (gSac_4factors locations)
% The 5th group (if present) is gSac_jph
n_locs = min(n_groups, 4);

%% Define feature names (19 features)
feature_names = { ...
    'best_vis', '2nd_vis', ...
    '3rd_vis', 'worst_vis', ...
    'best_sac', '2nd_sac', ...
    '3rd_sac', 'worst_sac', ...
    'best_delay', 'worst_delay', ...
    'n_loc_sig_vis', 'n_loc_sig_sac', ...
    'n_loc_sig_delay', 'n_loc_suppressed', ...
    'vis_to_sac_ratio', 'delay_index', ...
    'has_suppression', 'has_mixed_polarity', ...
    'tuning_breadth' ...
    };
nFeatures = length(feature_names);

%% Extract raw features for each neuron
raw_features = nan(nClusters, nFeatures);

for i_c = 1:nClusters

    %% Compute FR change from baseline at each loc
    % group_mean_frs{g}(cluster, epoch):
    %   epoch 1=baseline, 2=visual, 3=delay, 4=saccade
    fr_vis = nan(n_locs, 1);
    fr_sac = nan(n_locs, 1);
    fr_del = nan(n_locs, 1);

    for i_g = 1:n_locs
        bsl = group_mean_frs{i_g}(i_c, 1);
        fr_vis(i_g) = ...
            group_mean_frs{i_g}(i_c, 2) - bsl;
        fr_del(i_g) = ...
            group_mean_frs{i_g}(i_c, 3) - bsl;
        fr_sac(i_g) = ...
            group_mean_frs{i_g}(i_c, 4) - bsl;
    end

    %% Rank locations by visual response (descending)
    [~, rank_order] = sort(fr_vis, 'descend');

    fr_vis_ranked = fr_vis(rank_order);
    fr_sac_ranked = fr_sac(rank_order);
    fr_del_ranked = fr_del(rank_order);

    % Handle case where n_locs < 4
    % (pad with NaN if fewer than 4 locations)
    if n_locs < 4
        fr_vis_ranked(end+1:4) = NaN;
        fr_sac_ranked(end+1:4) = NaN;
        fr_del_ranked(end+1:4) = NaN;
    end

    %% Features 1-4: ranked visual FR changes
    raw_features(i_c, 1) = fr_vis_ranked(1);
    raw_features(i_c, 2) = fr_vis_ranked(2);
    raw_features(i_c, 3) = fr_vis_ranked(3);
    raw_features(i_c, 4) = fr_vis_ranked(4);

    %% Features 5-8: ranked saccade FR changes
    % (same rank order as visual)
    raw_features(i_c, 5) = fr_sac_ranked(1);
    raw_features(i_c, 6) = fr_sac_ranked(2);
    raw_features(i_c, 7) = fr_sac_ranked(3);
    raw_features(i_c, 8) = fr_sac_ranked(4);

    %% Features 9-10: best/worst delay FR changes
    raw_features(i_c, 9) = fr_del_ranked(1);
    raw_features(i_c, 10) = fr_del_ranked(4);

    %% Features 11-14: significance counts
    % sig_by_loc: [nClusters x nGroups x 3]
    %   dim3: 1=vis, 2=delay, 3=sac
    n_sig_vis = 0;
    n_sig_sac = 0;
    n_sig_del = 0;
    n_supp = 0;
    for i_g = 1:n_locs
        if sig_by_loc(i_c, i_g, 1)
            n_sig_vis = n_sig_vis + 1;
        end
        if sig_by_loc(i_c, i_g, 3)
            n_sig_sac = n_sig_sac + 1;
        end
        if sig_by_loc(i_c, i_g, 2)
            n_sig_del = n_sig_del + 1;
        end
        if any(supp_by_loc(i_c, i_g, :))
            n_supp = n_supp + 1;
        end
    end
    raw_features(i_c, 11) = n_sig_vis;
    raw_features(i_c, 12) = n_sig_sac;
    raw_features(i_c, 13) = n_sig_del;
    raw_features(i_c, 14) = n_supp;

    %% Features 15-16: temporal ratios
    eps_val = 0.1;
    best_vis_abs = abs(fr_vis_ranked(1));
    best_sac_abs = abs(fr_sac_ranked(1));
    best_del_abs = abs(fr_del_ranked(1));

    % vis_to_sac_ratio
    raw_features(i_c, 15) = best_vis_abs / ...
        (best_vis_abs + best_sac_abs + eps_val);

    % delay_index
    raw_features(i_c, 16) = best_del_abs / ...
        (max(best_vis_abs, best_sac_abs) + eps_val);

    %% Features 17-18: polarity
    any_supp = any( ...
        supp_by_loc(i_c, 1:n_locs, :), 'all');
    any_act = any( ...
        sig_by_loc(i_c, 1:n_locs, :), 'all');
    raw_features(i_c, 17) = double(any_supp);
    raw_features(i_c, 18) = ...
        double(any_act && any_supp);

    %% Feature 19: tuning breadth
    % best_vis / mean(all 4 vis changes)
    % High = sharply tuned, near 1 = omnidirectional
    mean_vis = mean(abs(fr_vis_ranked(1:4)), ...
        'omitnan');
    if mean_vis > eps_val
        raw_features(i_c, 19) = ...
            best_vis_abs / mean_vis;
    else
        % No visual response at any location
        raw_features(i_c, 19) = 1;
    end
end

%% Determine quality gate status
passes_quality = false(nClusters, 1);
for i_c = 1:nClusters
    is_modulated = ...
        any(sig_by_loc(i_c, 1:n_locs, :), ...
        'all') || ...
        any(supp_by_loc(i_c, 1:n_locs, :), ...
        'all');

    % Check valid FR data at any location
    has_valid_frs = false;
    for i_g = 1:n_locs
        if ~any(isnan( ...
                group_mean_frs{i_g}(i_c, :)))
            has_valid_frs = true;
            break;
        end
    end

    passes_quality(i_c) = ...
        is_modulated && has_valid_frs;
end

%% Normalize features (z-score across included)
feature_matrix = raw_features;
included_idx = find(passes_quality);

if ~isempty(included_idx)
    for i_f = 1:nFeatures
        vals = raw_features(included_idx, i_f);
        mu = mean(vals, 'omitnan');
        sd = std(vals, 'omitnan');
        if sd > 0
            feature_matrix(included_idx, i_f) = ...
                (vals - mu) / sd;
        else
            feature_matrix(included_idx, i_f) = 0;
        end
    end
end

%% Pack output
features = struct();
features.cluster_ids = cluster_ids;
features.feature_matrix = feature_matrix;
features.raw_features = raw_features;
features.feature_names = feature_names;
features.passes_quality = passes_quality;
features.scSide = scSide;

fprintf(['extract_sc_response_features: ' ...
    '%d/%d neurons pass quality gates. ' ...
    '%d features extracted.\n'], ...
    sum(passes_quality), nClusters, nFeatures);

end
