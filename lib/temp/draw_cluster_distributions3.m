function draw_cluster_distributions3(obj, C, output_prefix, labelsDimensions, listOfDimensionsToPrint, colorsForLines)
    numClusters = max(C)-1; %The last cluster is the mixed one, don't show here.
    numDimensionsToPrint = max(size(listOfDimensionsToPrint));
    
    columnsPlot = min(2, numDimensionsToPrint);
    rowsPlot = min(12, ceil(numDimensionsToPrint/2));

    for i=1:numDimensionsToPrint
        subplot(rowsPlot, columnsPlot, mod(i-1, 24)+1);
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
        if (mod(i-1, 24) == 23)
            fig = gcf;
            fig.PaperPosition = [0 0 20 20];
            fig_name = 'cluster_distribution_';
            print([output_prefix fig_name num2str(ceil(i/24)) '.tif'], '-dtiff');
            figure
        end
    end
end
