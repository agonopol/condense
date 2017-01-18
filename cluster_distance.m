function [distances, dist_significance] = cluster_distance(M, cluster_idx, dimensions, varargin)
    num_partitions = 20;
    compute_significance = false;
    num_samples = 100; 
    distfun = 'mi';
    num_datapoints = size(M, 1);
    
    for i=1:length(varargin)
        if(strcmp(varargin{i},'num_partitions'))
            num_partitions = varargin{i+1};
        end
        if(strcmp(varargin{i},'compute_significance'))
            compute_significance = varargin{i+1};
        end
        if(strcmp(varargin{i},'num_samples'))
            num_samples = varargin{i+1};
        end
        if(strcmp(varargin{i},'distfun'))
            distfun = varargin{i+1};
        end
    end
    
    num_dimensions = size(dimensions, 2);
    gene_ids =  zeros(1, num_dimensions);
    distances = zeros(1, num_dimensions);
    dist_significance = zeros(1, num_dimensions); 
    
    for i=1:num_dimensions
        if(strcmp(distfun,'L1'))
            data = M(:, dimensions(i));
            distances(i)=compute_L1_dist(data, cluster_idx, num_partitions, varargin{:});
        elseif(strcmp(distfun,'EMD'))
            data = M(:, dimensions(i));
            distances(i) = fast_emd_distro_1d(data, data(cluster_idx), varargin);
        elseif(strcmp(distfun, 'shift'))
            data = M(:, dimensions(i));
            distances(i) = compute_shift_dist(data, cluster_idx, varargin{:});
        else
            data = M(:, dimensions(i));
            distances(i)=compute_hist_differential_entropy(data, cluster_idx, num_partitions);
        end
    end
    
    if (compute_significance)
        for i=1:num_dimensions
            sample_distances = zeros(num_samples,1);
            if (strcmp(distfun,'L1'))
                for j=1:num_samples
                    sample_idx = randsample(num_datapoints,length(cluster_idx));
                    data = M(:, dimensions(i));
                    sample_distances(j)=compute_L1_dist(data, sample_idx, num_partitions);
                end
                dist_significance(i) = length(find(sample_distances>=distances(i)))/num_samples;
            elseif(strcmp(distfun,'EMD'))
                for j=1:num_samples
                    sample_idx = randsample(num_datapoints,length(cluster_idx));  
                    data = M(:, dimensions(i));
                    sample_distances(j) = fast_emd_distro_1d(data, data(sample_idx), varargin);
                end
                dist_significance(i) = length(find(sample_distances>= distances(i)))/num_samples;
            elseif(strcmp(distfun,'shift'))
                dataX = M(:, dimensions(i));
                dataY = M(cluster_idx, dimensions(i));
                if distances(i) > 0
                    dist_significance(i) = ranksum(dataX, dataY, 'tail', 'left');
                else
                    dist_significance(i) = ranksum(dataX, dataY, 'tail', 'right');
                end
            else
                for j=1:num_samples
                    sample_idx = randsample(num_datapoints,length(cluster_idx));
                    data = M(:, dimensions(i));
                    sample_distances(j)=compute_hist_differential_entropy(data, sample_idx, num_partitions);
                end
                dist_significance(i) = length(find(sample_distances >= distances(i)))/num_samples;
            end
        end
    end
end
