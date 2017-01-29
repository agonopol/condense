function drawClusterDistributions(samples, ...
                                  clusterAssignment, ...
                                  varargin)
% DRAWCLUSTERDISTRIBUTIONS Emits a plot with density estimates for dimensions
% of the data.
%
% Author : Tobias Welp, modifying a version presumably by Smita Krishnaswamy
% Year   : 2016

    % default arguments
    clustersToOmit = [];
    listOfDimensionsToEmit = 1:size(samples, 2);
    columnsPerPage = 1;
    rowsPerPage = 8;
    prefixFilename = '';

    % parsing of variable argument list
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'clustersToOmit'))
            clustersToOmit = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'colorsClusters'))
            colorsClusters = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'listOfDimensionsToEmit'))
            listOfDimensionsToEmit = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'labelsDimensions'))
            labelsDimensions = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'columnsPerPage'))
            columnsPerPage = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'rowsPerPage'))
            rowsPerPage = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'prefixFilename'))
            prefixFilename = varargin{i+1};
        end
    end

    % process arguments
    numClusters = max(clusterAssignment);
    assert(numClusters == length(unique(clusterAssignment)), ...
           'No proper format of cluster assignment');
    clustersToEmit = setdiff(1:numClusters, clustersToOmit);
    numClusters = length(unique(clustersToEmit));

    if (~exist('colorsClusters', 'var'))
        colorsClusters = distinguishable_colors(numClusters);
    end
    assert(length(colorsClusters) >= numClusters, ...
           'Insufficient number of colors passed');
    if (size(colorsClusters, 1) > numClusters)
        warning('More colors passed than required. Only using a subset of colors');
    end

    numDimensionsToEmit = length(listOfDimensionsToEmit);

    if (~exist('labelsDimensions'))
        labelsDimensions = {};
        for i=1:size(samples,2)
            if (any(i==listOfDimensionsToEmit))
                labelsDimensions = [labelsDimensions ; ['Dimension_' num2str(i)] ];
            end
        end
    end
    assert(length(labelsDimensions) >= length(listOfDimensionsToEmit), ...
           'Insufficient number of colors passed');
    if (length(labelsDimensions) > length(listOfDimensionsToEmit))
        warning('More dimension labels passed than required. Only using a subset of labels');
    end
    
    % start plotting cluster distributions.
    for i=1:numDimensionsToEmit
        subplot(rowsPerPage, columnsPerPage, i);
        hold on;
        dimensionData = samples(:, listOfDimensionsToEmit(i));
        [f,xi] = ksdensity(dimensionData);
        area(xi, f, 'FaceColor', [0.65,0.65,0.65], 'linestyle', 'none');
        maxDensity = max(f);
        for j=1:numClusters
            clusterIdx = (clusterAssignment==clustersToEmit(j));
            c = ksdensity(samples(clusterIdx, listOfDimensionsToEmit(i)), xi); 
            plot(xi, c, 'linewidth', 2, 'Color', colorsClusters(j, :));
            maxDensity = max([maxDensity c]);
        end
        xlabel(labelsDimensions(i));
        xlim([min(xi) max(xi)]);
        rangeY = [0 maxDensity*1.1];
        ylim(rangeY);
        line([mean(dimensionData) mean(dimensionData)], rangeY, 'Color', 'black', 'LineStyle', ':');
        for j=1:numClusters
            clusterIdx = (clusterAssignment==clustersToEmit(j));
            clusterMean = mean(samples(clusterIdx, listOfDimensionsToEmit(i)));
            line([clusterMean clusterMean], rangeY, 'Color', colorsClusters(j, :), 'LineStyle', '--');
        end
    end
    fig = gcf;
    fig.PaperPosition = [0 0 20 20];
    fig_name = '_cluster_distribution';
    print([prefixFilename fig_name '.tif'], '-dtiff');
end
