function [conflatedClusterAssignment] = conflateClusterAssignment(clusterAssignment, sampleIndices)
    conflatedClusterAssignment = [];
    for i = 1:length(sampleIndices)
        currentSampleIndices = sampleIndices{i};
        %%for j = 1:(length(currentSampleIndices)-1)
        %%    assert(isequal(clusterAssignment(currentSampleIndices(j)), ...
        %%                   clusterAssignment(currentSampleIndices(j+1))));
        %%end
        conflatedClusterAssignment(i, :) = clusterAssignment(currentSampleIndices(1), :);
    end
    assert(size(conflatedClusterAssignment, 1)==length(sampleIndices), ...
           'The number of samples in the conflated space must match that of the cluster assignment');
end
