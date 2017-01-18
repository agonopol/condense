function C = cluster_louvain(M, varargin)
k = 10;
npca = [];
for i=1:length(varargin)-1
    if(strcmp(varargin{i},'k'))
        k = varargin{i+1};
    end
    if(strcmp(varargin{i},'npca'))
        npca = varargin{i+1};
    end
end
if ~isempty(npca)
    disp 'PCA before distances and louvain'
    M = pcaMaaten(M, npca);
end
[D,I] = calcDistanceMatrix(M, 'k_knn', k, 'distfun', 'euclidean');
[A] = calcAffinityMatrix(D, I, 'sigma', Inf);
C = louvain(full(A));
end
