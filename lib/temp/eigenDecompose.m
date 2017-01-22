function [ eigenvectors, eigenvalues ] = eigenDecompose(M, varargin)
% EIGENDECOMPOSE Returns the eigendecomposition of the matrix supplied
%   per argument.
%
%   Author: Tobias Welp
%   Year  : 2016

    % default arguments
    mode = 'normal';
    nystroemN = -1;

    % parsing of variable argument list
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'mode'))
            mode = varargin{i+1};
            if ~ismember(mode, {'normal', 'nystroem'})
                warning(['eigenDecompose: Invalid choice for argument' ...
                         'mode. Using default value "normal"']);
                mode = 'normal';
            end
        end
        if (strcmp(varargin{i}, 'nystroemN'))
            nystroemN = varargin{i+1};
        end
    end

    if ((nystroemN ~= -1) & (mode ~= 'nystroem'))
        warning(['Set Nystroem parameter has no effect as nystroem ' ...
                 'has not been selected as mode']);
    end

    switch (mode)
        case 'normal'
            [eigenvectors, eigenvalues] = eigs(M);
        case 'nystroem'
            n = min(nystroemN, min(size(M)));
            A = full(M(1:n, 1:n));
            B = full(M(1:n, n+1:end));
            A_negSqrt = sqrtm(pinv(A));
            S = A + A_negSqrt * B * B' * A_negSqrt;
            [U, eigenvalues, T] = svd(S);
            eigenvectors = [A ; B'] * A_negSqrt * U * pinv(sqrt(eigenvalues));
    end
end
