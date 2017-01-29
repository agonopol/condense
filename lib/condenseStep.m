
function condensed = condenseStep(data)
    [distanceMatrix, indicatorMatrix] = calcDistanceMatrix(data);
    [normalizedAffinityMatrix, Q] = calcNormalizedAffinityMatrix(distanceMatrix, indicatorMatrix);
    diffusedNormalizedAffinityMatrix = diffuse(normalizedAffinityMatrix, Q);
    condensed =  diffusedNormalizedAffinityMatrix * data;
end