%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run_ODE_Pop_SMMALA.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% User Interface for simplified manifold MALA Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Run_ODE_Pop_SMMALA()

    clear all;  % Clear workspace 
    tic;        % Start timer

    % Script to add paths and do other things to enable algorithm to run
    startup();

    % Specify MCMC Options
    Options.burnin          = 2000;
    Options.thin            = 100   ;
    Options.numSamples      = 2000;
    Options.resultsName     = 'Model1';
    Options.numTemp         = 10;
    Options.invTemp         = cell(Options.numTemp, 1);
    tempVec                 = zeros(Options.numTemp);
    for i = 2:Options.numTemp
        tempVec(i) = 1 - 0.8^(Options.numTemp - i + 1);
        Options.invTemp{i} = tempVec(i);
    end
    Options.invTemp{1} = 1.0;
    
    Options.MaxACLag        = 20;
    Options.RandEps         = true;
    Options.MutateProb      = 0.66;

    % Function Handles for Model
    Functions.LogPosterior = @LogPostModel1;
    Functions.GibbsSampler = @GibbsSampleSigmaModel1;
    Functions.MetricTensor = @MetricTensorModel1;
    Functions.CreateParmStruct = @CreateParmStructModel1;

    % Data
    load('./data/model1_data.mat');
    stateData = [Data.A, Data.B];
    timePoints = t;

    % Initialize parameters
    numTemp = numel(Options.invTemp);
    w = 0.5;    % 0 for totally random, 1 for parm in simulation
    parameters = cell(numTemp, 1);
    for i = 1:numTemp
       parameters{i} = RandInitModel1(0); 
    end
    parmIndex = CreateParmIndex(parameters{1});
    randParmIndex = CreateRandParmIndex(parameters{1});

    % Specify MALA parameters   
    SMMALAparameters.StepSize     = 1e-2;    % Mean step size for Euler-Maruyama Integrator
    SMMALAparameters.ErrorBound   = 1;     % Upper limit on number of solver errors  

    % Run ODE_HMC
    ODE_Pop_SMMALA(stateData, timePoints, parameters, parmIndex, randParmIndex, ...
                Functions, Options, SMMALAparameters);

    % Stop timer
    toc;

end