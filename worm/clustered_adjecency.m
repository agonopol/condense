function [D, Z] = clustered_adjecency(data, clusters)
    D = zeros(length(unique(clusters)), length(unique(clusters)));
    Z = sparse(diag(ones(1, length(unique(clusters)))));
    for cluster = unique(clusters)
        for source = find(clusters == cluster)
            for sink = find(data(source, :))
                D(cluster, clusters(sink)) = D(cluster, clusters(sink)) + data(source, sink);
            end
        end
    end
    
end