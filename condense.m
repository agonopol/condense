close all;
clear;
clc;

files = dir('data/*.fcs');
for file = files'    
     hasDebris(fullfile(file.folder, file.name), fullfile(file.folder, strrep(file.name, 'splorm_0_normalized_Clean.fcs', 'DNA1_vs_DN2.png')));
end
