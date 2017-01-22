function regularizedClusterAssignment = regularizeClusterAssignment(clusterAssignment)
    regularizedClusterAssignment = clusterAssignment;
    for i = 1:length(clusterAssignment)
        while (sum(regularizedClusterAssignment == i) == 0)
            indicesToDecrement = find(regularizedClusterAssignment > i);
            if (isempty(indicesToDecrement))
                break;
            end
            regularizedClusterAssignment(indicesToDecrement) = ...
                regularizedClusterAssignment(indicesToDecrement)-1;
        end
    end
end
