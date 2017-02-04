classdef OptionsContractionClustering
    properties
        verbosityLevel = 1;
        indentationLevel = 0;
        numPrincipalComponents = 50;
        numDiffusionSteps = 1;
        maxNumContractionSteps = 60;
        kKNN = 20;
        initialSigma = Inf;
        nystroemN = 200;
        plotAnimation = true;
        inertia = 0.5;
        expertLabels = [];
        expId = '';
        emitRuntimeBreakdown = false;
        emitCondensationSequence = false;
        mergeEpsilonClusters = true;
        thresholdControlSigma = 0.01;
        clusterAssignmentMethod;
        controlSigmaMethod;
        frequencyMergingEpsilonClusters;
        modeCalcDistances = 'normal';
        epsilonClusterIdentificationMethod = 'constantEpsilon';
        prefixFileNames = '';
    end
    methods
        function obj = OptionsContractionClustering(varargin)
            mode = 'classic';
            for i=1:length(varargin)-1
                if (strcmp(varargin{i}, 'mode'))
                    mode = varargin{i+1};
                end
            end
            
            switch (mode)
                case 'classic'
                    obj.clusterAssignmentMethod = 'spectral';
                    obj.controlSigmaMethod = 'nuclearNormStabilization';
                    obj.frequencyMergingEpsilonClusters = 'uponMetastability';
                    
                case 'fast'
                    obj.clusterAssignmentMethod = 'none';
                    obj.controlSigmaMethod = 'movementStabilization';
                    obj.frequencyMergingEpsilonClusters = 'always';
            end
        end
        function str = asString(obj)
            str = ['cc_' ...
                   obj.expId '_' ...
                   num2str(obj.numPrincipalComponents) '_' ...
                   num2str(obj.numDiffusionSteps) '_' ...
                   num2str(obj.kKNN) '_' ...
                   num2str(obj.initialSigma) '_' ...
                   num2str(obj.nystroemN) '_' ...
                   num2str(obj.inertia) '_' ...
                   num2str(obj.mergeEpsilonClusters) '_' ...
                   obj.clusterAssignmentMethod '_' ...
                   obj.controlSigmaMethod '_' ...
                   obj.frequencyMergingEpsilonClusters '_' ...
                   obj.epsilonClusterIdentificationMethod '_' ...
                   obj.modeCalcDistances ...
                  ];
        end
    end
end
