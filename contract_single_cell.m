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
options.destination = fullfile(pwd(), 'results', 'contract', name, 'DB2_umi50_norm', '//');
[dest, ~, ~] = fileparts(options.destination);
mkdir_if_not_exists(dest);
 
%remove empty rows / columsn 
data.DB2_umi50_norm( ~any(data.DB2_umi50_norm,2), : ) = [];  %rows
data.DB2_umi50_norm( :, ~any(data.DB2_umi50_norm,1) ) = [];  %columns

contractor = ContractionClustering(data.DB2_umi50_norm, [], options);
contractor = contractor.contract();
clc;
