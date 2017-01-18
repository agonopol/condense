function [ resultM, resultSampleIndices ] = conflateClusters(M, sampleIndices, clusterAssignment)
    resultM = [];
    resultSampleIndices = {};
    numClusters = max(clusterAssignment);
    for i = 1:numClusters
        rowIndicesInCluster = find(clusterAssignment==i);
        clusterMedian = median(M(rowIndicesInCluster, :), 1);
        resultM = [resultM; clusterMedian];
        resultSampleIndices{i} = sampleIndices{rowIndicesInCluster(1)};
        for j = 2:length(rowIndicesInCluster)
            resultSampleIndices{i} = union(resultSampleIndices{i}, ...
                                           sampleIndices{rowIndicesInCluster(j)});
        end
    end
end
