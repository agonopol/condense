close all;
clear;
clc;

addpath('./lib');

% generate affinity matrix, noramlized
% diffuse

% load the data, clean transform
path = fullfile(pwd(), 'data/6300_Blood_D1_splorm_0_normalized_Clean.fcs');
obj = CyTOFData(path);
obj.dataTransformed = CyTOFData.transform(obj.data, 1);
markers = obj.markerNames';
channels = find(arrayfun(@(x1) ischar(x1{1}) && ~isempty(strfind(x1{1}, '_')) && isempty(strfind(x1{1}, 'DNA')) , markers));
data = obj.dataTransformed(:, channels');

% generate distance matrix
k = 10;
[distanceMatrix, indicatorMatrix] = calcDistanceMatrix(data, 'k_knn', k, 'distfun', 'euclidean');
[normalizedAffinityMatrix, Q] = calcNormalizedAffinityMatrix(distanceMatrix, indicatorMatrix, 'sigma', Inf);
diffusedNormalizedAffinityMatrix = diffuse(normalizedAffinityMatrix, Q);
% C = louvain(full(diffusedNormalizedAffinityMatrix));