function draw_cluster_distributions2(obj, C, cluster_sig_genes, cluster_sig_signs, max_genes, output_prefix, cluster_sigs, labelsDimensions, listOfDimensionsToPrint, colorsForLines)
    numClusters = max(C)-1; %The last cluster is the mixed one, don't show here.
    numDimensionsToPrint = max(size(listOfDimensionsToPrint));
    
    columnsPlot = min(2, numDimensionsToPrint);
    rowsPlot = min(6, ceil(numDimensionsToPrint/2));

    for i=1:numDimensionsToPrint
        subplot(rowsPlot, columnsPlot, mod(i-1, 12)+1);
        hold on;
        dimensionData = obj.data(:, listOfDimensionsToPrint(i));
        [f,xi] = ksdensity(dimensionData);
        area(xi,f,'FaceColor',[0.65,0.65,0.65], 'linestyle', 'none');
        maxDensity = max(f);
        for j=1:numClusters
            clusterIdx = (C==j);
            c = ksdensity(obj.data(clusterIdx, listOfDimensionsToPrint(i)), xi); 
            plot(xi, c, 'linewidth', 2, 'Color', colorsForLines(j, :));
            maxDensity = max([maxDensity c]);
        end
        xlabel(labelsDimensions(listOfDimensionsToPrint(i)));
        xlim([min(xi) max(xi)]);
        rangeY = [0 maxDensity*1.1];
        ylim(rangeY);
        line([mean(dimensionData) mean(dimensionData)], rangeY,'Color', 'black', 'LineStyle', '--');
        for j=1:numClusters
            clusterIdx = (C==j);
            clusterMean = mean(obj.data(clusterIdx, listOfDimensionsToPrint(i)));
            line([clusterMean clusterMean], rangeY, 'Color', colorsForLines(j, :), 'LineStyle', '--');
        end
        if (mod(i-1, 12) == 11)
            fig_name = '_cluster_distribution_';
            print([output_prefix fig_name num2str(ceil(i/12))], '-dtiff');
            figure
        end
    end

    %%[rows, cols] = size(cluster_sig_genes);

    %%for i=1:num_clusters
    %%    figure;
    %%    if(cols == 1)
    %%        current_cluster_genes = cluster_sig_genes{i}
    %%        current_cluster_signs = cluster_sig_signs{i}
    %%    else
    %%        current_cluster_genes = cluster_sig_genes(i,:)
    %%        current_cluster_signs = cluster_sig_signs(i,:)
    %%    end
    %%    
    %%    cluster_idx = (C==i);
    %%    ngenes = min(max_genes, length(current_cluster_genes));
    %%    for j=1:ngenes
    %%        nr = ceil(sqrt(ngenes));
    %%        nc = nr;
    %%        subplot(nr, nc, j);
    %%        hold on;
    %%        current_gene = current_cluster_genes(j);
    %%        current_sign = current_cluster_signs(j);
    %%        cluster_data = obj.data(cluster_idx,current_gene);
    %%        [f,xi] = ksdensity(obj.data(:,current_gene));
    %%        
    %%        area(xi,f,'FaceColor',[0.65,0.65,0.65], 'linestyle', 'none');
    %%        %set(gca,'FontSize',14);
    %%        
    %%        c = ksdensity(cluster_data, xi); 
    %%        plot(xi,c,'linewidth', 2, 'Color','r');
    %%        if(current_sign>0)
    %%            xlabel(strcat(labelsDimensions(current_gene),'+'));
    %%        elseif(current_sign<0)
    %%            xlabel(strcat(labelsDimensions(current_gene),'-'));
    %%        else
    %%            xlabel(labelsDimensions(current_gene));
    %%        end
    %%        xlim([min(xi) max(xi)]);
    %%        ylim([0 max([c f])*1.1]);
    %%        %ylabel('density');
    %%        %draw lines at the means
    %%        cluster_mean = mean(cluster_data);
    %%        population_mean = mean(obj.data(:,current_gene));
    %%        YL = ylim;
    %%        line([cluster_mean cluster_mean],[0 YL(2)],'Color','r');
    %%        line([population_mean population_mean],[0 YL(2)],'Color', 'k');
    %%        if exist('cluster_sigs','var')
    %%            title(['p = ' num2str(cluster_sigs{i}(j),2)]);
    %%        end
    %%    end
    %%    
    %%    titlestring = sprintf('cluster %d',i);
    %%    mtit(titlestring);
    %%    
    %%    nr = ceil(sqrt(ngenes));
    %%    nc = nr;
    %%    fig_name = 'cluster_distribution_';
    %%    set(gcf, 'PaperPosition', [0, 0, nc * 3, nr * 2]);
    %%    print('-dtiff',[output_prefix '_' fig_name num2str(i) '.tiff']);
    %%    
    %%end
end
