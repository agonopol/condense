function [] = hasDebris(path, out)

obj = CyTOFData(path);
markers = obj.markerNames';
dna1 = obj.data(:, arrayfun(@(x1) ischar(x1{1}) && ~isempty(strfind(x1{1}, 'DNA1')) , markers));
dna2 = obj.data(:, arrayfun(@(x1) ischar(x1{1}) && ~isempty(strfind(x1{1}, 'DNA2')) , markers));

scatter(asinh(dna1/5), asinh(dna2/5));
print(out, '-dpng');
end