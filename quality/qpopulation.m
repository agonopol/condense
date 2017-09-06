function qpopulation(X, Y, limit, measure, varargin)
   output = 'quality.png';
   
   if (strcmp(measure, 'mmd'))
        quality = @mmdscore;
   elseif (strcmp(measure, 'emd'))
        quality = @emdscore;
   elseif (strcmp(measure, 'rindex'))
        quality = @rindexscore;
   end
   
   neurons = intersect(Y.channels, X.channels);
   iX = getindex(X.channels, neurons);
   iY = getindex(Y.channels, neurons);
   
   k = intersect(arrayfun(@(x) max(X.clusterAssignments(x, :)), 1:size(X.clusterAssignments, 1)), ...
                arrayfun(@(y) max(Y.clusterAssignments(y, :)), 1:size(Y.clusterAssignments, 1)));
   k  = k(k < limit);
   pX = arrayfun(@(i) any(k == i), arrayfun(@(x) max(X.clusterAssignments(x, :)), 1:size(X.clusterAssignments, 1)));
   pX = orderP(unique(X.clusterAssignments(pX, :), 'rows'));
   
   pY = arrayfun(@(i) any(k == i), arrayfun(@(y) max(Y.clusterAssignments(y, :)), 1:size(Y.clusterAssignments, 1)));
   pY = orderP(unique(Y.clusterAssignments(pY, :), 'rows'));
   
   dX = X.contractionSequence(:, :, 1);
   dY = Y.contractionSequence(:, :, 1);
   
   hold on;
   ylabel(measure);

   [x, y] = quality(dX, dY, pX, pY, iX, iY);
   plot(x, y, 'DisplayName' , 'condensation');
   
   for i=1:2:length(varargin)-1
        if (strcmp(varargin{i}, 'output'))
            output = varargin{i+1};
        else
            label = varargin{i};
            fn = varargin{i+1};
            [x, y] = quality(dX, dY, fn(dX, pX), fn(dY, pY), iX, iY);
            plot (x, y, 'DisplayName' , label);
        end
   end
   
   xlabel('# of clusters') % x-axis label
   title('PairWise Clustering Quality')

   legend ( 'show' );
   hold off;

   saveas(gcf, output);
   close all;
   
end

function [x, scores] = mmdscore(dX, dY, pX, pY, iX, iY)
    x = arrayfun(@(x) max(pX(x, :)), 1:size(pX, 1));
    scores = zeros(size(x));
    for i = 1:size(pX, 1)
       [cx, ~] = stats(dX, pX(i, :)); 
       [cy, ~] = stats(dY, pY(i, :));
       distance = mmd(cx, cy);
       if isreal(distance)
            scores(i) = abs(distance);
       else
            scores(i) = nan;
        end
    end
end

function [P] = orderP(P)
    [~, ii] = sort(max(P'));
    P = P(ii, :);
end

function [x, scores] = rindexscore(dX, dY, pX, pY, iX, iY)
    x = arrayfun(@(x) max(pX(x,:)), 1:size(pX, 1));
    scores = zeros(size(x));
    for i = 1:size(pX, 1)
       distance = rindex(pX(i, iX), pY(i, iY));
       scores(i) = abs(distance);
    end
end

function [x, scores] = emdscore(dX, dY, pX, pY, iX, iY)
    x = arrayfun(@(x) max(pX(x,:)), 1:size(pX, 1));
    scores = zeros(size(x));
    for i = 1:size(pX, 1)
       [cx, wx] = stats(dX, pX(i, :)); 
       [cy, wy] = stats(dY, pY(i, :));
       distance = emd(cx, cy, (wx ./ sum(wx))', (wy ./ sum(wy))');
       if isreal(distance)
            scores(i) = abs(distance);
       else
            scores(i) = nan;
        end
    end
end

function index = getindex(intersection, neurons)
    index = arrayfun(@(y) ...
        find(arrayfun(@(x) ~isempty(x{:}), strfind(intersection, y))), neurons, 'UniformOutput', false);
    index = cell2mat(index);
    index = index';
end