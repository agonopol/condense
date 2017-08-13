function q = affmodularity(data, assigment)    
    % calculate affinity matrix
    A = CalculateAffinity(data);
    q = modularity(A, assigment);
end

function [affinity] = CalculateAffinity(data)
    % set the parameters
    sigma = 1;
    D = squareform(pdist(data));

    for i=1:size(data,1)    
        for j=1:size(data,1)
            dist = D(i,j);
            affinity(i,j) = exp(-dist/(2*sigma^2));
        end
    end
end