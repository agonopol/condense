function [I, J] = kneepoint(X, Y, varargin)
    fn = @mmd_mesh;
    maxclusters = 1000;
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'distance'))
            if (strcmp(varargin{i+1}, 'mmd'))
                fn = @mmd_mesh;
            else
                fn = @emd_mesh;
            end
        end
        if (strcmp(varargin{i}, 'maxclusters'))
            maxclusters = varargin{i+1};
        end
    end
    A = area(X, Y, maxclusters);
    [D] = fn(X,Y, A);
    mesh(D);
    [I, J] = knee(D);
end

function [I, J] = knee(D)
    [I, J] = find(diff(D) > 0, 1, 'first');
end

function [M] = area(X, Y, num)
    i = find( arrayfun(@(x) length(unique(X.clusterAssignments(x, :))), 1:size(X.clusterAssignments, 1)) < num, 1, 'first');
    n = find( arrayfun(@(x) length(unique(X.clusterAssignments(x, :))), 1:size(X.clusterAssignments, 1)) < 2, 1, 'first');
    
    j = find( arrayfun(@(x) length(unique(Y.clusterAssignments(x, :))), 1:size(Y.clusterAssignments, 1)) < num, 1, 'first');
    m = find( arrayfun(@(x) length(unique(Y.clusterAssignments(x, :))), 1:size(Y.clusterAssignments, 1)) < 2, 1, 'first');
    
    M = NaN(n, m);
    for x = i:n
        for y = j:m
            M(x, y) = 0.0;
        end
    end
end

function [MMD] = mmd_mesh(X, Y, MMD)
    dX = X.contractionSequence(:, :, 1);
    dY = Y.contractionSequence(:, :, 1);

    for x = 1:size(MMD, 1)
        [cx, wx] = stats(dX, X.clusterAssignments(x, :));
        for y = 1:size(MMD, 2)
            if ~isnan(MMD(x,y))
                [cy, ~] = stats(dY, Y.clusterAssignments(y, :));
                distance = mmd(cx, cy);
                if isreal(distance)
                    MMD(x,y) = abs(distance);
                else
                    MMD(x, y) = nan;
                end
            end
        end
    end
end

function [EMD] = emd_mesh(X, Y, EMD)

    dX = X.contractionSequence(:, :, 1);
    dY = Y.contractionSequence(:, :, 1);

    for x = 1:size(EMD, 1)
        [cx, wx] = stats(dX, X.clusterAssignments(x, :));
        for y = 1:size(EMD, 2)
            if ~isnan(EMD(x,y))
                [cy, wy] = stats(dY, Y.clusterAssignments(y, :));
                distance = emd(cx, cy, (wx ./ sum(wx))', (wy ./ sum(wy))');
                if isreal(distance)
                    EMD(x, y) = abs(distance);
                else
                    EMD(x, y) = nan;
                end
            end
        end
    end
end