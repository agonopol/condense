function combined = combineStep( data )
    [distanceMatrix, ~] = calcDistanceMatrix(data);
    epsilon = mean2(distanceMatrix) - std2(distanceMatrix);
    [i,j, ~] = find((distanceMatrix < epsilon) & tril(ones(size(distanceMatrix)), -1));
    combined = data;
    for index=1:length(i)
        row = mean([data(i(index), :); data(j(index), :)]);
        combined(i(index), :) = row;
    end
    combined(setdiff(j,i), :) = [];
end

