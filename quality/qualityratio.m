function [ ratios ] = qualityratio( scores, assigments )
    ratios = zeros(size(assigments));
    for i = 2:size(assigments, 1)
        for j = 1:max(assigments(i, :))
            current = score(scores, assigments, i, j);
            previous = 0.0;
            pclusters = unique(assigments(i - 1, assigments(i, :) == j));
            if size(pclusters, 2) == 1
                ratios(i, assigments(i, :) == j) = ratios(i - 1, find(assigments(i, :) == j, 1, 'first'));
            else
                for k = 1:size(pclusters, 2)
                    pc = pclusters(k);
                    previous = previous + score(scores, assigments, i - 1, pc);
                end
                if previous > 0.0 
                    ratios(i, assigments(i, :) == j) = current / previous;
                else
                    ratios(i, assigments(i, :) == j) = 0.0;
                end
            end
        end
    end
end

function q = score( scores, assigments, level, cluster)
    index = find(assigments(level, :) == cluster, 1, 'first');
    q = scores(level, index);
end