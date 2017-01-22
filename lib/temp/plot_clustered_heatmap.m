function plot_clustered_heatmap(obj, C, labelsDimensions, varargin)
% plot_clustered_heatmap(obj, C, varargin)

M = obj.data;
% This removes the last cluster. Do not want to do this any longer.
%mixedCluster = max(C); 
%indicesMixedCluster = C==mixedCluster;
%M = M(find(~(C==mixedCluster)),:);
%C = C(find(~(C==mixedCluster)));

clim = [];
channel_names = labelsDimensions;
if(isprop(obj,'genes'))
    channel_names = obj.genes;
end
if(isprop(obj,'channel_name_map'))
    channel_names = obj.channel_name_map;
end

cmp = parula(64);
cluster_rows = false;
cluster_cols = false;
cluster_rows_mean = false;
cluster_cols_mean = false;
cluster_rows_sort = false;
%cmp_clusters = lines(size(unique(C),2)); %%'lines';
cmp_clusters = distinguishable_colors(size(unique(C),2)); %%'lines';
sort_vec = [];
sort_vec_no_bar = [];
col_first = [];
columnClusters = [];
zscore_cols = false;
clusterMeanPerSample = false;
for i=1:length(varargin)-1
    if(strcmp(varargin{i},'data'))
        M = varargin{i+1};
    end
    if (strcmp(varargin{i},'clusterMeanPerSample'))
        clusterMeanPerSample=varargin{i+1};
    end
    if(strcmp(varargin{i},'colormap'))
        cmp = varargin{i+1};
    end
    if(strcmp(varargin{i},'sort_vec'))
        sort_vec = varargin{i+1};
    end
    if(strcmp(varargin{i},'sort_vec_no_bar'))
        sort_vec_no_bar = varargin{i+1};
    end
    if(strcmp(varargin{i},'colormapclusters'))
        cmp_clusters = varargin{i+1};
    end
    if(strcmp(varargin{i},'columnClusters'))
        columnClusters = varargin{i+1};
    end
    if(strcmp(varargin{i},'clusterrows'))
        cluster_rows = varargin{i+1};
    end
    if(strcmp(varargin{i},'clustercols'))
        cluster_cols = varargin{i+1};
    end
    if(strcmp(varargin{i},'clusterrowsmean'))
        cluster_rows_mean = varargin{i+1};
    end
    if(strcmp(varargin{i},'clustercolsmean'))
        cluster_cols_mean = varargin{i+1};
    end
    if(strcmp(varargin{i},'clusterrowssort'))
        cluster_rows_sort = varargin{i+1};
    end
    if(strcmp(varargin{i},'clim'))
        clim = varargin{i+1};
    end
    if(strcmp(varargin{i},'col_first'))
        col_first = varargin{i+1};
    end
    if(strcmp(varargin{i},'zscore_cols'))
        zscore_cols = varargin{i+1};
    end
    if(strcmp(varargin{i}, 'channels'))
        sprintf('channels specified \n')
        [~, col_ind] = ismember(varargin{i+1}, channel_names);
        M = obj.data(:,col_ind(col_ind > 0));
    end
    if(strcmp(varargin{i}, 'channel_labels'))
        sprintf('channel names specified \n')
        channel_ind = ismember(varargin{i+1}, channel_names);
        channel_names = varargin{i+1};
        channel_names = channel_names(channel_ind);
    end
end

if (clusterMeanPerSample)
    % This stuff averages the expressions on each dimensions across all samples belonging
    % to a cluster.
    for i=1:size(unique(C),2)
        M(find(C==i),:) = ones(size(M(find(C==i),:)))*diag(mean(M(find(C==i),:), 1));
    end
end

C = C(:);
if ~isempty(sort_vec_no_bar)
    [~,ind] = sort(sort_vec_no_bar);
    M = M(ind,:);
end
if ~isempty(sort_vec)
    [sort_vec,ind] = sort(sort_vec);
    M = M(ind,:);
end
[C,ind] = sort(C);
M = M(ind,:);
if ~isempty(sort_vec)
    sort_vec = sort_vec(ind);
end
% cluster M
if cluster_rows_sort
    disp 'cluster rows sort'
    D_rows = squareform(pdist(M));
    Z_rows = linkage(D_rows);
    ind_rows = optimalleaforder(Z_rows,D_rows);
    M = M(ind_rows,:);
    C = C(ind_rows);
    [C,ind] = sort(C);
    M = M(ind,:);
end
if cluster_cols
    D_cols = squareform(pdist(M'));
    Z_cols = linkage(D_cols);
    ind_cols = optimalleaforder(Z_cols,D_cols);
    M = M(:,ind_cols);
    labelsDimensions = labelsDimensions(ind_cols);
end
if cluster_rows
    disp 'cluster rows'
    D_rows = squareform(pdist(M));
    Z_rows = linkage(D_rows);
    size(D_rows)
    size(Z_rows)
    ind_rows = optimalleaforder(Z_rows,D_rows);
    M = M(ind_rows,:);
    C = C(ind_rows);
end
if cluster_cols_mean
    disp 'cluster cols mean'
    M_mean = nan(size(M));
    for I=unique(C)'
        ind = C == I;
        M_mean(ind,:) = repmat(mean(M(ind,:)),sum(ind),1);
    end
    D_cols = squareform(pdist(M_mean'));
    Z_cols = linkage(D_cols);
    size(D_cols)
    size(Z_cols)
    ind_cols = optimalleaforder(Z_cols,D_cols);
    M = M(:,ind_cols);
    channel_names = channel_names(ind_cols);
end
if cluster_rows_mean
    disp 'cluster rows mean'
    M_mean = nan(size(M));
    for I=unique(C)'
        ind = C == I;
        M_mean(ind,:) = repmat(mean(M(ind,:)),sum(ind),1);
    end
    D_rows = squareform(pdist(M_mean));
    Z_rows = linkage(D_rows);
    size(D_rows)
    size(Z_rows)
    ind_rows = optimalleaforder(Z_rows,D_rows);
    M = M(ind_rows,:);
    C = C(ind_rows);
end
if ~isempty(col_first)
    disp 'force first col'
    ind_first = find(ismember(channel_names, col_first));
    M = M(:,[ind_first 1:ind_first-1 ind_first+1:end]);
    channel_names = channel_names([ind_first 1:ind_first-1 ind_first+1:end]);
end
if zscore_cols
    M = zscore(M);
end
for q = 1:size(M,2)
    M(:,q) = mat2gray(M(:,q));
end
figure;
% bar 1
ax1 = subplot(1,40,1);
colormap(ax1, cmp_clusters);
%colormap(lines(4));
imagesc(C);
hold on;
% white lines
Y_white = find(abs(diff(C)) > 0);
plot(repmat(xlim, length(Y_white), 1)', [Y_white(:) Y_white(:)]' + 0.5, '-w', 'linewidth', 2);
set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
%set(gca,'Visible','off')
% bar 2
if ~isempty(sort_vec)
    ax2 = subplot(1,40,2);
    colormap(ax2, cmp_clusters);
    imagesc(sort_vec);
    hold on;
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    % black lines
    Y_black = find(abs(diff(sort_vec)) > 0);
    plot(repmat(xlim, length(Y_black), 1)', [Y_black(:) Y_black(:)]' + 0.5, '-k', 'linewidth', .5);
    % white lines
    Y_white = find(abs(diff(C)) > 0);
    plot(repmat(xlim, length(Y_white), 1)', [Y_white(:) Y_white(:)]' + 0.5, '-w', 'linewidth', 2);
end
% main
if ~isempty(sort_vec)
    ax3 = subplot(1,40,3:40);
else
    ax3 = subplot(1,40,2:40);
end
colormap(ax3, cmp);
if ~isempty(clim)
    imagesc(M, clim);
else
    imagesc(M);
end
hold on;
if ~isempty(sort_vec)
    % black lines
    Y_black = find(abs(diff(sort_vec)) > 0);
    plot(repmat(xlim, length(Y_black), 1)', [Y_black(:) Y_black(:)]' + 0.5, '-k', 'linewidth', .5);
end
% white lines
Y_white = find(abs(diff(C)) > 0);
plot(repmat(xlim, length(Y_white), 1)', [Y_white(:) Y_white(:)]' + 0.5, '-w', 'linewidth', 2);
X_white = find(abs(diff(columnClusters)) > 0);
plot([X_white(:) X_white(:)]' + 0.5, repmat(ylim, length(X_white), 1)', '-w', 'linewidth', 2);
set(gca,'xtick',1:length(labelsDimensions));
set(gca,'xticklabel', labelsDimensions);
set(gca,'XTickLabelRotation',45);
set(gca,'TickLength',[ 0 0 ]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
xlhand = get(gca, 'xlabel');
set(gca, 'FontSize', 6);
fig = gcf;
fig.PaperPosition = [0 0 20 10];
colorbar;
%xticklabel_rotate([],45,[]);
pos_main = get(gca,'position');
pos_bar1 = get(ax1,'position');
pos_bar1([2 4]) = pos_main([2 4]);
set(ax1,'position',pos_bar1);
if ~isempty(sort_vec)
    pos_bar2 = get(ax2,'position');
    pos_bar2([2 4]) = pos_main([2 4]);
    set(ax2,'position',pos_bar2);
end
end
