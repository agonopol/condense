function [q, scores] = affmodularity(data, assigment)    
    % calculate affinity matrix
    [d, I] = calcDistanceMatrix(data);
    A = calcAffinityMatrix(d, I);
    
    w = sum(sum(A));
    scores = arrayfun(@(c) Q(A, assigment, w, c), 1:max(assigment));
    q = (1/(2*w)) * sum(scores);
end

function [q] = Q(A, assigment, w, cluster)
    q = 0;
    
    for i = find(assigment == cluster)
        wi = sum(A(i, :));
        for j = find(assigment == cluster)
            if ( i == j )
                continue
            end
            wj = sum(A(j, :));
            wij = A(i, j);
            q = q + (wij - ((wi * wj) / ( 2 * w)));
        end
    end
end