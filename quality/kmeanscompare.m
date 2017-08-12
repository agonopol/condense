function [assigments] = kmeanscompare(contracted)
    matrix = contracted.contractionSequence(:, :, 1);
    assigments = arrayfun(@(k) kmeans(matrix, max(contracted.clusterAssignments(k, :))), ...
        1:size(contracted.clusterAssignments, 1), ...
        'UniformOutput', false);
    assigments = cell2mat(assigments)';
end
