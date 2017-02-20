close all;
clear;
clc;

addpath('./lib');

% load the data, clean transform
path = fullfile(pwd(), 'data/6170_D1_splorm_0_normalized_Clean.fcs');
obj = CyTOFData(path);
obj.dataTransformed = CyTOFData.transform(obj.data, 1);
markers = obj.markerNames';
channels = find(arrayfun(@(x1) ischar(x1{1}) && ~isempty(strfind(x1{1}, '_')) && isempty(strfind(x1{1}, 'DNA')) , markers));
data = obj.dataTransformed(:, channels');
data = datasample(data, min(length(data), 1600),'Replace', false);

options = OptionsContractionClustering();
options.clusterAssignmentMethod = 'none';
options.controlSigmaMethod = 'always';
[~, name, ~] = fileparts(path);
options.destination = fullfile(pwd(), 'results', 'contract', name, '//');
[dest, ~, ~] = fileparts(options.destination);
mkdir_if_not_exists(dest);

contractor = ContractionClustering(data, options);
contractor = contractor.contract();

