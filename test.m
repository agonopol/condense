%% Main function to generate tests
function tests = test
    tests = functiontests(localfunctions);
end

%% Test Functions
function TestRetina(testCase)
   [data, channels] = retina('./test/dec24_sdata_raw.mat');
      
    options.destination = './test/output/retina//';
    
    [dest, ~, ~] = fileparts(options.destination);
    mkdir_if_not_exists(dest);
    
    contractor = ContractionClustering(data, channels, options);
    contractor.contract();
end


%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
    addpath('./test');
end


%% Optional fresh fixtures  
function setup(testCase)  % do not change function name

    options = Options();
    options.clusterAssignmentMethod = 'none';
    options.frequencyMergingEpsilonClusters = 'uponMetastability'; %always,uponMetastability%
    options.controlSigmaMethod = 'nuclearNormStabilization';
    options.numDiffusionSteps = 3;
    
end

function teardown(testCase)  % do not change function name
    
    clc;
    close all force;
    close all hidden;
    
end