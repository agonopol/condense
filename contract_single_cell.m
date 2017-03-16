close all;
clear;
clc;

addpath('./lib');

options = OptionsContractionClustering();
options.clusterAssignmentMethod = 'none';
options.frequencyMergingEpsilonClusters = 'uponMetastability'; %always,uponMetastability%
options.controlSigmaMethod = 'nuclearNormStabilization';
options.numDiffusionSteps = 1;
options.fastStop = false;

path = 'library_size_normalized_scRNA-seq_DB.mat';
data = load(path);
    
[~, name, ~] = fileparts(path);
options.destination = fullfile(pwd(), 'results', 'contract', name, 'DB1_umi50_norm', '//');
[dest, ~, ~] = fileparts(options.destination);
mkdir_if_not_exists(dest);
 
contractor = ContractionClustering(data.DB1_umi50_norm, [], options);
contractor = contractor.contract();
clc;
