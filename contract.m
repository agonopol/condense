close all;
clear;
clc;

addpath('./lib');

options = OptionsContractionClustering();
options.clusterAssignmentMethod = 'none';
options.frequencyMergingEpsilonClusters = 'uponMetastability';
options.controlSigmaMethod = 'movementStabilization';
options.numDiffusionSteps = 1;

files = dir('data/*.fcs');
for file = files'
    close all;
    clc;
    
    path = fullfile(file.folder, file.name);
    obj = CyTOFData(path);
    obj.dataTransformed = CyTOFData.transform(obj.data, 1);
    markers = obj.markerNames';
    channels = find(arrayfun(@(x1) ischar(x1{1}) && ~isempty(strfind(x1{1}, '_')) && isempty(strfind(x1{1}, 'DNA')) , markers));
    data = obj.dataTransformed(:, channels');
    data = datasample(data, min(length(data), 1600),'Replace', false);
    
    [~, name, ~] = fileparts(path);
    options.destination = fullfile(pwd(), 'results', 'contract', name, '//');
    [dest, ~, ~] = fileparts(options.destination);
    mkdir_if_not_exists(dest);
 
    contractor = ContractionClustering(data, options);
    contractor = contractor.contract();
end

close all;
clc;
