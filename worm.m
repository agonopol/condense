close all;
clear;
clc;

addpath('./lib');
addpath('./PHATE/Matlab');
addpath('./worm');
addpath('./lib/sanky/template');
addpath('./lib/sanky');


options = OptionsContractionClustering();
options.clusterAssignmentMethod = 'none';
options.epsilonClusterIdentificationMethod = 'constantEpsilon';
options.frequencyMergingEpsilonClusters = 'always'; %always,uponMetastability%
options.controlSigmaMethod = 'nuclearNormStabilization'; %nuclearNormStabilization,movementStabilization
options.numDiffusionSteps = 3;
options.fastStop = true;
options.maxClusters = 7;
options.phateEmbedding = false;

files = dir('worm/*.csv');


for file = files'
    path = fullfile(file.folder, file.name);
    data = readworm(path);
    [neurons, data] = weightedadj(data);
    
    [~, name, ~] = fileparts(path);
    options.destination = fullfile(pwd(), 'results', 'worm', name, '//');
    [dest, ~, ~] = fileparts(options.destination);
    mkdir_if_not_exists(dest);
    
    markov = rownorm(data);
    [V,D] = eig(markov);
    matrix = bsxfun(@times,V(:,1:50)',diag(D(1:50, :)))';
    
    contractor = ContractionClustering(matrix, cellstr(neurons), options);
    contractor = contractor.contract();
    
    
    f = fopen(strcat(options.asString(), '_assigments.csv'), 'w');
    
    fprintf(f, join(string(neurons'), ','));
    fprintf(f, '\n');
    flipped = flip(contractor.clusterAssignments);

    for i = 1:min(size(contractor.clusterAssignments))
        fprintf(f, join(string(flipped(i,:)), ','));
        fprintf(f, '\n');
    end
    fclose(f);

    sanky(contractor.clusterAssignments, neurons, strcat(options.asString(), '_sanky.html'));
    clc;
    close all force;
    close all hidden;
end