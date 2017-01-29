close all;
clear;
clc;

addpath('./lib');

% generate affinity matrix, noramlized
% diffuse

% load the data, clean transform
path = fullfile(pwd(), 'data/6300_HF_D1_splorm_0_normalized_Clean.fcs');
obj = CyTOFData(path);
obj.dataTransformed = CyTOFData.transform(obj.data, 1);
markers = obj.markerNames';
channels = find(arrayfun(@(x1) ischar(x1{1}) && ~isempty(strfind(x1{1}, '_')) && isempty(strfind(x1{1}, 'DNA')) , markers));
data = obj.dataTransformed(:, channels');

v = VideoWriter('results/condense/6300_HF_D1_tsne.avi');
v.FrameRate = 15;
open(v);

scatterX(data, 'dimensionalityReductionMethod', 'tsne');
% xlim([-12,12]);
% ylim([-12,12]);
writeVideo(v,getframe(gcf));
    
N=5;
for i = 1:N
    data = condenseStep(data);
    scatterX(data, 'dimensionalityReductionMethod', 'tsne');
%     xlim([-12,12]);
%     ylim([-12,12]);
    writeVideo(v,getframe(gcf));
end

close(v);
close all;
