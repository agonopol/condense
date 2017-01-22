function thresholdHistogram(data, threshold, varargin)

    colorAssignment = [0,0,1; 1,0,0];
    plotFractions = false;
    verticalBars = true;
    numBins = false;

    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'colorAssignment'))
            colorAssignment = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'plotFractions'))
            plotFractions = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'verticalBars'))
            verticalBars = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'numBins'))
            numBins = varargin{i+1};
        end
    end

    if (logical(numBins))
        h = histogram(data, numBins);
    else
        h = histogram(data);
    end
    values = h.Values;
    if (plotFractions)
        values = values/sum(values);
    end
    binEdges = h.BinEdges(2:end);
    valuesHigh = values .* (binEdges > threshold);
    valuesLow  = values .* (binEdges <= threshold);
    if (verticalBars)
        b = bar(binEdges, valuesHigh, 'FaceColor', colorAssignment(1,:), 'LineStyle', 'none');
    else
        b = barh(binEdges, valuesHigh, 'FaceColor', colorAssignment(1,:), 'LineStyle', 'none');
    end
    hold on;
    if (verticalBars)
        b = bar(binEdges, valuesLow, 'FaceColor', colorAssignment(2,:), 'LineStyle', 'none');
    else
        b = barh(binEdges, valuesLow, 'FaceColor', colorAssignment(2,:), 'LineStyle', 'none');
    end
    hold off;
end
