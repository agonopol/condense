function q = modularity(adj, s)

deg = sum(adj);
m = sum(deg)/2;
s = s(:);

q = s' * modmat(adj) * s / (4*m);

end

function [b] = modmat(adj)
    n = length(adj);
    deg = sum(adj,1);
    m = sum(deg)/2;
    b = zeros(size(adj));
    for k1 = 1:n
        for k2 = 1:n
            b(k1, k2) = adj(k1, k2) - deg(k1) * deg(k2) / (2 * m);
        end
    end
end