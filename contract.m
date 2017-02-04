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

options = OptionsContractionClustering();
contractor = ContractionClustering(data, options);
contractor.contract();

