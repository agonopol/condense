close all;
clear;
clc;

addpath('./lib');

options = OptionsContractionClustering();
options.clusterAssignmentMethod = 'none';
options.frequencyMergingEpsilonClusters = 'uponMetastability';
options.controlSigmaMethod = 'nuclearNormStabilization';
options.numDiffusionSteps = 3;
options.fastStop = true;

files = dir('data/*.fcs');
for file = files'
    close all;
    clc;
    
    path = fullfile(file.folder, file.name);
    obj = CyTOFData(path);
    obj.dataTransformed = CyTOFData.transform(obj.data, 1);
    fields = channels(obj);
    data = obj.dataTransformed(:, cell2mat(fields(:,1))');
    data = datasample(data, min(length(data), 1600),'Replace', false);
    
    [~, name, ~] = fileparts(path);
    options.destination = fullfile(pwd(), 'results', 'contract', name, '//');
    [dest, ~, ~] = fileparts(options.destination);
    mkdir_if_not_exists(dest);
 
    contractor = ContractionClustering(data, fields(:,2), options);
    contractor = contractor.contract();
    close all;
    contractor.plotClusterHeatMap();
end

close all;
clc;
