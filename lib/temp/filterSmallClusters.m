function [ filteredDataPoints, filteredClusterAssignment ] = filterSmallClusters(dataPoints, clusterAssignment, cutoff)
    % remove clusters smaller than cutoff data points
    clusters = unique(clusterAssignment);
    for i=1:length(clusters)
        if (sum(clusterAssignment == clusters(i)) < cutoff)
            dataPoints = dataPoints(clusterAssignment ~= clusters(i), :);
            clusterAssignment = clusterAssignment(clusterAssignment ~= clusters(i));
        end
    end

    % regularlize (cluster assignments are from 1 to #clusters)
    remainingClusters = unique(clusterAssignment);
    for i=1:length(remainingClusters)
        clusterAssignment(clusterAssignment == remainingClusters(i)) = i;
    end

    filteredDataPoints = dataPoints;
    filteredClusterAssignment = clusterAssignment;
end
