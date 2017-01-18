function [inflatedM] = inflateClusters(M, sampleIndices)
    inflatedM = [];
    for i = 1:size(M, 1)
        currentSampleIndices = sampleIndices{i};
        inflatedM(currentSampleIndices, :) = repmat(M(i, :),length(currentSampleIndices),1);
    end
end
