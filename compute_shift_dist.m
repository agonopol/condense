function mean_distance = compute_shift_dist(data, cluster_idx, varargin)
    mean_distance = mean(data(cluster_idx)) - mean(data);
end
