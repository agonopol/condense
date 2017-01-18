function [sig_cluster_dimensions, cluster_gene_signs, sig_dimensions, dimensions, cluster_sigs, gene_sigs] = cluster_discriminating_dimensions(M, clusterAssignment, varargin)
        
    sig_value = 5;
    sig_threshold = 0.05;
    num_clusters = size(unique(clusterAssignment), 1);
    compute_significance = false;
    gene_sigs = [];
    cluster_sigs = [];
    
    dimensions = [];
    sig_dimensions=[];
    
    for i=1:length(varargin)
        if(strcmp(varargin{i},'compute_significance'))
            compute_significance =  varargin{i+1};
        end
        if(strcmp(varargin{i},'sig_value'))
            sig_value =  varargin{i+1};
        end
        if(strcmp(varargin{i},'sig_threshold'))
            sig_threshold =  varargin{i+1};
        end
        if(strcmp(varargin{i},'dimensions'))
            dimensions =  varargin{i+1};
        end
    end
    
    %everything higher than a significance value
    sig_cluster_dimensions = cell(num_clusters,1);
    cluster_gene_signs = cell(num_clusters,1);
    cluster_sigs = cell(num_clusters,1);
    
    [distances, dist_significance, shifts, shift_significance, dimensions] = ...
        population_cluster_distance(M, ...
                                    clusterAssignment, ...
                                    varargin{:}, ...
                                    'compute_significance', true);

    for i=1:num_clusters
        sig_clusters_dimensions_idx = find(dist_significance(i,:)<sig_threshold);
       
        sig_cluster_dimensions{i} =  dimensions(sig_clusters_dimensions_idx);
        cluster_gene_signs{i} = shifts(i, sig_clusters_dimensions_idx);
        cluster_sigs{i} = dist_significance(i, sig_clusters_dimensions_idx);
        
        %returning in order of most significant
        [cluster_sigs{i}, sorted_idx] = sort(cluster_sigs{i});
        sig_cluster_dimensions{i} = sig_cluster_dimensions{i}(sorted_idx);
        cluster_gene_signs{i} = cluster_gene_signs{i}(sorted_idx); 
    end
     
    gene_avg_sig = mean(dist_significance);
    sig_gene_ids = gene_avg_sig<sig_threshold;
    sig_dimensions = dimensions(sig_gene_ids);
    gene_sigs = gene_avg_sig(sig_gene_ids);
    
    %returning in order of most significant
    [gene_sigs, sorted_idx] = sort(gene_sigs);
    sig_dimensions = sig_dimensions(sorted_idx);
    %going to return this in sorted order
end
