function [distances, dist_significance, shifts, shift_significance, dimensions ] = population_cluster_distance(M, clusterAssignment, varargin)
        
    distfun = 'mi';
    compute_significance = false;
    num_clusters = size(unique(clusterAssignment), 1);
    dimensions = [];
    
    for i=1:length(varargin)
       if(strcmp(varargin{i},'dimensions'))
           dimensions = varargin{i+1};
       end
    end
    num_dimensions = size(dimensions, 2);
    
    distances = zeros(num_clusters, num_dimensions);
    dist_significance = zeros(num_clusters, num_dimensions);
    shifts = zeros(num_clusters, num_dimensions);
    shift_significance = zeros(num_clusters, num_dimensions);
    
    for i=1:length(varargin)
        if(strcmp(varargin{i},'distfun'))
            distfun = varargin{i+1};
        end
        if(strcmp(varargin{i},'compute_significance'))
            compute_significance = varargin{i+1};
        end
    end
    
    for i=1:num_clusters
        cluster_idx = find(clusterAssignment == i);
        [distances(i,:), dist_significance(i,:)] = cluster_distance(M, cluster_idx, dimensions, varargin{:});
    end

    for i=1:num_clusters
        cluster_idx = find(clusterAssignment == i);
        
        [shifts(i,:), shift_significance(i,:)] = cluster_distance(M, cluster_idx, dimensions, 'distfun', 'shift','compute_significance', true);
       
        posshifts = (shifts(i,:)>0)&(shift_significance(i,:)<.05);
        negshifts = (shifts(i,:)<0)&(shift_significance(i,:)<.05);
        negshifts = negshifts*(-1);
        
        shifts(i,:) = posshifts+negshifts;
    end
    
end
