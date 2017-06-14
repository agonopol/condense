function affinityMatrix = calcAffinityMatrix(distanceMatrix, zeroIndicatorMatrix, varargin)
% CALCAFFINITYMATRIX Returns affinity matrix corresponding to
%   the distance matrix given as argument.
%
%   Author: Tobias Welp, modifying a version by Ingo Buerk
%   Year  : 2016 based on original work from 2011/2012

    % default arguments
    sigma = 1;
    exponent = 2;

    % parsing of variable argument list
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'sigma'))
            sigma = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'exponent'))
            exponent = varargin{i+1};
        end
    end

    affinityMatrix = spfun(@(distanceMatrix) (simGaussian(distanceMatrix, sigma, exponent)), distanceMatrix);
    affinityMatrix = affinityMatrix + zeroIndicatorMatrix;
end

function [ M ] = simGaussian( M, sigma, k )
%SIMGAUSSIAN Calculates Gaussian similarity on matrix
%   simGaussian(M, sigma) returns a matrix of the same size as
%   the distance matrix M, which contains similarity values
%   that are computed by using a Gaussian similarity function
%   with parameter sigma.
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

    if nargin < 3
        k = 2;
    end

    M = exp(-M.^k ./ (2*sigma^2));

end
