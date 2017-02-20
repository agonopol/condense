classdef ContractionClustering
    properties
        options = [];
        % This holds the data point positions at each point of the sequence.
        contractionSequence = [];
        % This holds the data point positions at the current contraction step with
        % epsilon clusters being represented by a single point been merged.
        dataContracted = [];
        eigenvalueSequence = [];
        currentEigenvectors = [];
        currentEigenvalues = [];
        clusterAssignments = [];
        currentSigma = Inf;
        sampleIndices = [];
        normalizedAffinityMatrix = [];
        normalizedAffinityMatrixInitialSigma = [];
        weights = [];
        runtimes;
        iteration = 1;
        epsilon;
    end
    methods
        function obj = ContractionClustering(data, options)
            obj.options = options;
            obj.dataContracted = data;
            obj.currentSigma = obj.options.initialSigma;
            if (obj.currentSigma == Inf)
                [distanceMatrix, ~] = calcDistanceMatrix(data);
                obj.currentSigma = mean2(distanceMatrix) - std2(distanceMatrix);
            end

            % Make sure the expert labels are properly populated.
            if (isempty(obj.options.expertLabels))
                obj.options.expertLabels = ones(size(obj.dataContracted, 1), 1);
            end
            assert(isequal(length(obj.options.expertLabels), size(obj.dataContracted, 1)));

            % Populate sample indices.
            obj.sampleIndices = num2cell(1:size(obj.dataContracted, 1));

            % Calculate epsilon
            obj.epsilon = max(max(obj.dataContracted)-min(obj.dataContracted))/10000;

            obj.runtimes = containers.Map({'aff', 'spd', 'clas', 'vis', 'contr', 'rest'}, zeros(1, 6));
        end
        function obj = contract(obj)
            for iteration = 1:obj.options.maxNumContractionSteps
                obj.iteration = iteration;
                obj.contractionSequence(:, :, obj.iteration) = inflateClusters(obj.dataContracted, obj.sampleIndices);
                obj = obj.calcAffinities();
                obj = obj.spectralDecompose();
                obj = obj.performContractionMove();
                obj = obj.mergeEpsilonClusters();
                obj = obj.assignClusters();
                obj = obj.controlSigma();
                obj = obj.plotAlgorithmState();
                obj.printProgress(false);
                if (obj.checkTerminationCondition())
                    break;
                end
            end
            obj.emitRuntimeBreakdown();
            obj.emitCondensationSequence();
        end
        function obj = calcAffinities(obj)
            tic;
            [D, Z] = calcDistanceMatrix(obj.dataContracted, ...
                                        'k_knn', obj.options.kKNN, ...
                                        'type_k_knn', 'normal', ...
                                        'distfun', 'euclidean', ...
                                        'lengthPartitions', 2*obj.currentSigma, ...
                                        'mode', obj.options.modeCalcDistances, ...
                                        'n_pca', obj.options.numPrincipalComponents, ...
                                        'indentationLevel', obj.options.indentationLevel + 1, ...
                                        'verbosityLevel', obj.options.verbosityLevel - 1);
            obj.weights = cellfun(@length, obj.sampleIndices);
            if (obj.requireSpectralDecomposition())
                obj.normalizedAffinityMatrixInitialSigma = calcNormalizedAffinityMatrix(D, Z, ...
                                                                                        'sigma', obj.options.initialSigma, ...
                                                                                        'exponent', 2, ...
                                                                                        'weights', obj.weights);
            end
            obj.normalizedAffinityMatrix = calcNormalizedAffinityMatrix(D, Z, ...
                                                                        'sigma', obj.currentSigma, ...
                                                                        'exponent', 2, ...
                                                                        'weights', obj.weights);
            obj.runtimes('aff') = obj.runtimes('aff') + toc;
        end
        function obj = spectralDecompose(obj)
            if (obj.requireSpectralDecomposition())
                tic
                [obj.currentEigenvectors, obj.currentEigenvalues] = ...
                    eig(obj.normalizedAffinityMatrixInitialSigma*diag(obj.weights)+0.00000001);
                [obj.currentEigenvalues, indicesSort] = sort(abs(diag(obj.currentEigenvalues)));
                obj.eigenvalueSequence(obj.iteration, :) = ...
                    [min(min(obj.eigenvalueSequence))*ones(size(obj.contractionSequence, 1)-length(obj.currentEigenvalues), 1) ; obj.currentEigenvalues];
                obj.currentEigenvectors = obj.currentEigenvectors(:, indicesSort);
                obj.currentEigenvectors = abs(obj.currentEigenvectors(:, obj.currentEigenvalues > 0.99));
                obj.runtimes('spd') = obj.runtimes('spd') + toc;
            end
        end
        function obj = assignClusters(obj)
            tic;
            switch (obj.options.clusterAssignmentMethod)
                case 'none'
                    clusterAssignment = (1:size(obj.dataContracted, 1))';
                case 'spectral'
                    %% Runs spectral clustering. Note that this means that k-means
                    %% is run on the eigenvectors corresponding to the eigenvalues
                    %% close to one. I doubt that it is technically spectral clustering
                    %% because this would require that we did take the eigenvectors
                    %% of the laplacian and we do not really use the laplacian here.
                    if (size(obj.currentEigenvectors, 2) > 20)
                        obj.currentEigenvectors = pcaMaaten(obj.currentEigenvectors, 20);
                    end
                    clusterAssignment = kmeans(obj.currentEigenvectors, sum(obj.currentEigenvalues>0.99));
                otherwise
                    error(['Unknown cluster assignment method: ' obj.options.clusterAssignmentMethod]);
            end
            obj.clusterAssignments(obj.iteration, :) = inflateClusters(clusterAssignment, ...
                                                                       obj.sampleIndices);
            obj.runtimes('clas') = obj.runtimes('clas') + toc;
        end
        function obj = plotAlgorithmState(obj)
            persistent myTextObj;
            persistent lastSigma;
            persistent sigmaBumps;
            persistent currentLengthIterations;
            persistent relativeMovement;
            persistent previousDataContracted;
            if (obj.options.plotAnimation)
                tic
                if obj.iteration == 1
                    sigmaBumps = [];
                    lastSigma = obj.currentSigma;
                    set(gcf, 'Position', [2068 1 1200 800]);
                    currentLengthIterations = 50;
                    relativeMovement = [];
                    previousDataContracted = obj.dataContracted;
                else
                    if (obj.currentSigma ~= lastSigma)
                        lastSigma = obj.currentSigma;
                        sigmaBumps = [sigmaBumps obj.iteration];
                    end
                    if (obj.iteration > currentLengthIterations)
                        currentLengthIterations = 2*currentLengthIterations;
                    end
                end
                % Plotting samples at original position.
                ax1 = subplot('Position', [0.05, 0.40, 0.425, 0.425]);
                scatterX(obj.contractionSequence(:, :, 1), 'colorAssignment', obj.clusterAssignments(obj.iteration, :));
                colormap(ax1, distinguishable_colors(length(unique(obj.clusterAssignments(obj.iteration, :)))));
                lims = axis;
                % Plotting samples at contracted position.
                ax2 = subplot('Position', [0.525, 0.40, 0.425, 0.425]);
                sizeAssignment = sqrt(cellfun(@size, obj.sampleIndices, repmat({2}, 1, length(obj.sampleIndices))));
                scatterX(obj.dataContracted, ...
                         'colorAssignment', conflateClusterAssignment(obj.options.expertLabels, obj.sampleIndices), ...
                         'sizeAssignment', sizeAssignment');
                colormap(ax2, distinguishable_colors(length(unique(obj.options.expertLabels))));
                ax3 = subplot('Position', [0.05, 0.10, 0.9, 0.225]);
                switch (obj.options.controlSigmaMethod)
                    case 'nuclearNormStabilization'
                        % Plotting eigenvalue map.
                        colormap(ax3, parula(64));
                        eigenvalueMap = repmat(min(min(obj.eigenvalueSequence)), size(obj.eigenvalueSequence, 2), currentLengthIterations);
                        eigenvalueMap(1:size(obj.eigenvalueSequence, 2), 1:size(obj.eigenvalueSequence, 1)) = obj.eigenvalueSequence';
                        imagesc(eigenvalueMap);
                        for j = 1:size(sigmaBumps, 2)
                            line([sigmaBumps(j) sigmaBumps(j)], ylim, 'Color', 'red', 'LineWidth', 2);
                        end
                    case 'movementStabilization'
                        if (obj.iteration == 1)
                            relativeMovement = [ NaN ];
                        elseif (isequal(size(obj.dataContracted), size(previousDataContracted)))
                            relativeMovement = [relativeMovement ...
                                                max(sum(abs(obj.dataContracted-previousDataContracted)))/max(max(obj.dataContracted)-min(obj.dataContracted))+eps];
                        else
                            relativeMovement = [relativeMovement NaN];
                        end
                        previousDataContracted = obj.dataContracted;
                        plot(relativeMovement, 'o');
                        set(gca, 'YScale', 'log')
                        ylim([eps, max(max(relativeMovement), 1)])
                        set(gca, 'YTick', sort([eps, get(gca, 'YTick')]));
                        yticks = get(gca, 'YTick');
                        ytickLabels = {'eps'};
                        for i = 2:length(yticks)
                            ytickLabels{i} = num2str(yticks(i));
                        end
                        set(gca,'yticklabel', ytickLabels);
                        line([0 currentLengthIterations], [obj.options.thresholdControlSigma obj.options.thresholdControlSigma], 'LineStyle', '--')
                        for j = 1:size(sigmaBumps, 2)
                            line([sigmaBumps(j) sigmaBumps(j)], ylim, 'Color', 'red', 'LineWidth', 2);
                        end
                end
                % Plotting Header Line
                subplot('Position', [0.05, 0.85, 0.9, 0.125], 'Visible', 'off')
                toWrite = ['Iteration ' num2str(obj.iteration) ...
                           ', \sigma = ' num2str(obj.currentSigma) ...
                           ', #Clusters = ' num2str(length(unique(obj.clusterAssignments(obj.iteration, :)))) ...
                           ', #Samples = ' num2str(size(obj.dataContracted, 1))];
                if obj.iteration == 1
                    myTextObj = text(0, 0.5, toWrite, 'FontSize', 20);
                else
                    set(myTextObj, 'String', toWrite);
                end
                obj.saveFigureAsAnimationFrame();
                obj.runtimes('vis') = obj.runtimes('vis') + toc;
            end
        end
        function rsl = checkTerminationCondition(obj)
            tic
            rsl = false;
            numClusters = length(unique(obj.clusterAssignments(obj.iteration, :)));
            if (numClusters == 1 && obj.iteration > 10)
                rsl = true;
            end
            obj.runtimes('rest') = obj.runtimes('rest') + toc;
        end
        function obj = performContractionMove(obj)
            tic
            diffusedNormalizedAffinityMatrix = diffuse(obj.normalizedAffinityMatrix, 'numDiffusionSteps', obj.options.numDiffusionSteps, 'weights', obj.weights); 
            obj.dataContracted =   (1-obj.options.inertia) * weightedMultiply(diffusedNormalizedAffinityMatrix, obj.dataContracted, obj.weights) ...
                                 +  obj.options.inertia    * obj.dataContracted;
            obj.runtimes('contr') = obj.runtimes('contr') + toc;
        end
        function obj = controlSigma(obj)
            tic
            persistent iterationLastIncrease;
            persistent previousDataContracted;
            if (obj.iteration == 1)
                iterationLastIncrease = 1;
                previousDataContracted = obj.dataContracted;
            else
                switch (obj.options.controlSigmaMethod)
                    case 'nuclearNormStabilization'
                        %% The idea of this mode is to keep the sigma constant until the sum of eigenvalues
                        %% stabilizes. The sum of eigenvalues is considered stablized if
                        %% the sum of ten consecutive decreases is less than 5% of the
                        %% total sum of eigenvalues. After the sum of eigenvalues stabilized,
                        %% the sigma is increased by 20%.
                        if (   (obj.iteration > 10) ...
                            && (  sum(sum(obj.eigenvalueSequence(:, obj.iteration-10:obj.iteration-1))-sum(obj.eigenvalueSequence(:, obj.iteration-9:obj.iteration))) ...
                                < 0.05 * sum(obj.eigenvalueSequence(:, obj.iteration-10))))
                            obj.currentSigma = 1.1*obj.currentSigma;
                            disp(['Bumped sigma in iteration ' num2str(obj.iteration)]);
                        end
                    case 'movementStabilization'
                        if (isequal(size(obj.dataContracted), size(previousDataContracted)))
                            thisRelativeMovement = max(sum(abs(obj.dataContracted-previousDataContracted)))/max(max(obj.dataContracted)-min(obj.dataContracted))+eps;
                            disp(['Did not contract, checking if should bump sigma on itration ' num2str(obj.iteration), ...
                                ' with relative movment of ', num2str(thisRelativeMovement), '<', num2str(obj.options.thresholdControlSigma) ]);
                            if (thisRelativeMovement < obj.options.thresholdControlSigma)
                                obj.currentSigma = 1.1*obj.currentSigma;
                                disp(['Bumped sigma to ', num2str(obj.currentSigma), 'in iteration ', num2str(obj.iteration), ... 
                                    ' previous bump was on ', num2str(iterationLastIncrease)]);
                                iterationLastIncrease = obj.iteration;
                            end
                        end
                        previousDataContracted = obj.dataContracted;
                end
            end
            obj.runtimes('rest') = obj.runtimes('rest') + toc;
        end
        function obj = mergeEpsilonClusters(obj)
            tic
            persistent previousSigma;
            if (obj.iteration == 1)
                previousSigma = obj.currentSigma;
            else
                mergeEpsilonClusters = false;
                switch (obj.options.frequencyMergingEpsilonClusters)
                    case 'uponMetastability'
                        mergeEpsilonClusters = (obj.currentSigma ~= previousSigma);
                        previousSigma = obj.currentSigma;
                    case 'always'
                        mergeEpsilonClusters = true;
                end
                if (mergeEpsilonClusters)
                    disp('Merging Clusters');
                    switch (obj.options.epsilonClusterIdentificationMethod)
                        case 'constantEpsilon'
                            epsilon = obj.epsilon;
                        case 'dynamicSigmaFraction'
                            epsilon = obj.currentSigma/4;
                    end
                    [obj.dataContracted, obj.sampleIndices] = ...
                        conflateClusters(obj.dataContracted, ...
                                         obj.sampleIndices, ...
                                         detectEpsilonClusters(obj.dataContracted, epsilon));
                end
            end
            obj.runtimes('rest') = obj.runtimes('rest') + toc;
        end
        function rsl = requireSpectralDecomposition(obj)
            rsl = (   strcmp(obj.options.clusterAssignmentMethod, 'spectral') ...
                   || strcmp(obj.options.controlSigmaMethod, 'nuclearNormStabilization'));
        end
        function printProgress(obj, forcePrint)
            persistent timeLastPrint;
            if (obj.iteration==1)
                timeLastPrint = 0;
            end
            overallTime = sum(cell2mat(obj.runtimes.values())) + eps;
            if (forcePrint || overallTime > timeLastPrint + 10)
                timeLastPrint = overallTime;
                if (obj.options.verbosityLevel > 0)
                    message = ['ContractionClustering: Iteration ' sprintf('%4u', obj.iteration) ...
                               ', runtime: ' sprintf('%7.2f', overallTime)];
                    for key = obj.runtimes.keys()
                        message = [message ' ' key{1} '=' num2str(obj.runtimes(key{1})/overallTime, '%.2f')];
                    end
                    disp([indent(obj.options.indentationLevel) message]);
                end
            end
        end
        function emitRuntimeBreakdown(obj)
            if (obj.options.emitRuntimeBreakdown)
                runtimeLabels = obj.runtimes.keys();
                runtimes = obj.runtimes.values(runtimeLabels);
                save([obj.options.prefixFileNames obj.options.asString() '_runtimeBreakdown.mat'], ...
                     'runtimes', 'runtimeLabels');
            end
        end
        function emitCondensationSequence(obj)
            if (obj.options.emitCondensationSequence)
                condensationSequence = obj.contractionSequence;
                save([obj.options.prefixFileNames obj.options.asString() '_condensationSequence.mat'], ...
                     'condensationSequence');
            end
        end
        function saveFigureAsAnimationFrame(obj)
            im = print('-RGBImage');
            [imind,cm] = rgb2ind(im,256);
            filename = strcat(obj.options.asString(), '_animation.gif');
            if obj.iteration == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
            end
        end
        function emitClusterResults(obj)
            dlmwrite(strcat(obj.options.prefixFileNames, obj.options.asString(), '_clusterAssignments.txt'), obj.clusterAssignments);
        end
    end
end
