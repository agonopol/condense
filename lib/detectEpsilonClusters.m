function [ epsilonClusterAssignment ] = detectEpsilonClusters(M, epsilon)
    numSamples = size(M, 1);
    idx = knnsearch_fast(M, M, size(M, 1)-1);
    epsilonClusterAssignment = 1:numSamples;
    numClusters=1;
    for i=1:numSamples
        if (epsilonClusterAssignment(i) == i)
            epsilonClusterAssignment(i) = numClusters;
            for j=idx(i, :)
                if (norm(M(i, :)-M(j, :))<epsilon)
                    epsilonClusterAssignment(j) = numClusters;
                else
                    break;
                end
            end
            numClusters = numClusters+1;
        end
    end
end
