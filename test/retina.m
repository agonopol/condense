%% Read DataSet
function [data, channels] = retina(file)

sdata_raw =  scRNA_data('file_name', file);
%% filter lib size
sdata_sub = sdata_raw;
ind_keep = sdata_raw.library_size >= 1000;
sum(ind_keep)
sdata_sub.data = sdata_sub.data(ind_keep,:);
sdata_sub.cells = sdata_sub.cells(ind_keep);
genes_keep = sum(sdata_sub.data) > 0;
sdata_sub.data = sdata_sub.data(:, genes_keep);
sdata_sub.genes = sdata_sub.genes(genes_keep);
sdata_sub.mpg = sdata_sub.mpg(genes_keep);
sdata_sub.cpg = sdata_sub.cpg(genes_keep);
sdata_sub = sdata_sub.recompute_name_channel_map();
exp_vec = exp_vec(ind_keep);

%% get experiment
%sdata_sub = sdata_raw;
ind_keep = strcmp(exp_vec, 'Bipolar5');
sum(ind_keep)
sdata_sub.data = sdata_sub.data(ind_keep,:);
sdata_sub.cells = sdata_sub.cells(ind_keep);
genes_keep = sum(sdata_sub.data) > 0;
sdata_sub.data = sdata_sub.data(:, genes_keep);
sdata_sub.genes = sdata_sub.genes(genes_keep);
sdata_sub.mpg = sdata_sub.mpg(genes_keep);
sdata_sub.cpg = sdata_sub.cpg(genes_keep);
sdata_sub = sdata_sub.recompute_name_channel_map()

%% lib size norm global
sdata = sdata_sub;
sdata.data_nn = sdata.data;
sdata = sdata.normalize_data_fix_zero();

%% log scale
sdata.data = log(sdata.data + 0.1);

data = sdata.data;
end