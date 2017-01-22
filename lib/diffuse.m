    function [diffusedNormalizedAffinityMatrix] = diffuse(normalizedAffinityMatrix, varargin)

    numDiffusionSteps = 50;
    mode = 'exact';
    nystroemN = 200;
    rsvd = 10;

    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'numDiffusionSteps'))
            numDiffusionSteps = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'mode'))
            mode = varargin{i+1};
            if ~ismember(mode, {'exact', 'nystroem', 'rsvd', 'svd', 'eig'})
                warning(['Diffuse: Invalid choice "' mode '" for argument ' ...
                         'mode. Using default value "exact"']);
                mode = 'exact';
            end
        end
        if (strcmp(varargin{i}, 'nystroemN'))
            nystroemN = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'rsvdK'))
            rsvdK = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'weights'))
            weights = varargin{i+1};
        end
    end

    normalizedAffinityMatrix = full(normalizedAffinityMatrix);
    
    assert(numDiffusionSteps >= 1);

    switch (mode)
        case 'exact'
            if (exist('weights', 'var'))
                y = diag(1./weights);
                while (numDiffusionSteps > 1)
                    if (mod(numDiffusionSteps, 2)) 
                        %odd
                        y = weightedMultiply(y, normalizedAffinityMatrix, weights);
                        normalizedAffinityMatrix = weightedMultiply(normalizedAffinityMatrix, normalizedAffinityMatrix, weights);
                        numDiffusionSteps = (numDiffusionSteps-1)/2;
                    else
                        %even
                        normalizedAffinityMatrix = weightedMultiply(normalizedAffinityMatrix, normalizedAffinityMatrix, weights);
                        numDiffusionSteps = (numDiffusionSteps-1)/2;
                    end
                end
                diffusedNormalizedAffinityMatrix = weightedMultiply(y, normalizedAffinityMatrix, weights);
            else
                diffusedNormalizedAffinityMatrix = normalizedAffinityMatrix^numDiffusionSteps;
            end
        case 'svd'
            assert(~exist('weights', 'var'));
            [U,S,V] = svd(normalizedAffinityMatrix);
            Sd = S^numDiffusionSteps;
            diffusedNormalizedAffinityMatrix = U*Sd*V';
        case 'eig'
            assert(~exist('weights', 'var'));
            [U,L] = eigenDecompose(normalizedAffinityMatrix, ...
                                   'mode', 'normal');
            Ld = L^numDiffusionSteps;
            diffusedNormalizedAffinityMatrix = U*Ld*inv(U);
        case 'rsvd'
            assert(~exist('weights', 'var'));
            [U,S,V] = rndsvd(normalizedAffinityMatrix, rsvdK);
            Sd = S^numDiffusionSteps;
            diffusedNormalizedAffinityMatrix = U*Sd*V';
        case 'nystroem'
            assert(~exist('weights', 'var'));
            [ eigenvectors, eigenvalues ] = eigenDecompose(normalizedAffinityMatrix, ...
                                                           'mode', 'nystroem', ...
                                                           'nystroemN', nystroemN);
            diffusedEigenvalues = eigenvalues^numDiffusionSteps;
            diffusedNormalizedAffinityMatrix = eigenvectors*diffusedEigenvalues*eigenvectors';
    end
end
