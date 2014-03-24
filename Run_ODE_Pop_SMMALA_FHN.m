%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run_ODE_Pop_SMMALA_FHN.m
%
% Author: Chris Burke
% Last Modified: 02-11-14
%
% User Interface for simplified manifold MALA Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Run_ODE_Pop_SMMALA_FHN()

    clear all;  % Clear workspace 
    tic;        % Start timer

    % Script to add paths and do other things to enable algorithm to run
    startup();

    % Specify MCMC Options
    Options.burnin          = 0;
    Options.thin            = 1;
    Options.numSamples      = 2000;
    Options.resultsName     = 'Model1';
    Options.numTemp         = 3;
    
    Options.invTemp         = cell(Options.numTemp, 1);
    Options.invTemp{1}      = 1.0;
    Options.invTemp{2}      = 0.5;
    Options.invTemp{3}      = 0.1;
    
    Options.MaxACLag        = 20;
    Options.RandEps         = true;
    Options.MutateProb      = 0.9;

    % Function Handles for Model
    Functions.LogPosterior = @LogPostFHN;
    Functions.GibbsSampler = @GibbsSampleSigmaFHN;
    Functions.MetricTensor = @MetricTensorFHN;
    Functions.CreateParmStruct = @CreateParmStructFHN;

    % Data
    load('./data/model1_data.mat');
    stateData = [Data.A, Data.B];
    timePoints = t;

    % Initialize parameters
    numTemp = numel(Options.invTemp);
    w = 0.9;
    parameters = cell(numTemp, 1);
    for i = 1:numTemp
       parameters{i} = RandInitFHN(w); 
    end
    parmIndex = CreateParmIndex(parameters{1});
    randParmIndex = CreateRandParmIndex(parameters{1});

    % Specify MALA parameters   
    SMMALAparameters.StepSize     = 1e0;    % Mean step size for Euler-Maruyama Integrator
    SMMALAparameters.ErrorBound   = 1;     % Upper limit on number of solver errors  

    % Run ODE_HMC
    ODE_Pop_SMMALA(stateData, timePoints, parameters, parmIndex, randParmIndex, ...
                Functions, Options, SMMALAparameters);

    % Stop timer
    toc;

end