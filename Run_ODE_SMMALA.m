%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run_ODE_SMMALA.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% User Interface for simplified manifold MALA Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Run_ODE_SMMALA()

    clear all;  % Clear workspace 
    tic;        % Start timer

    % Script to add paths and do other things to enable algorithm to run
    StartupODEInference();

    % Specify MCMC Options
    Options.burnin          = 0;
    Options.thin            = 1;
    Options.numSamples    = 100;
    Options.resultsName     = 'Model1';
    Options.MaxACLag        = 20;
    Options.Monitor         = 100;

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
    weight = 0.5;
    parameters = RandInitModel1(weight);
    parmIndex = CreateRandParmIndex(parameters);

    % Specify MALA parameters   
    SMMALAparameters.StepSize     = 1e-8;    % Mean step size for Euler-Maruyama Integrator
    SMMALAparameters.ErrorBound   = 0.8;     % Upper limit on number of solver errors  

    % Run ODE_HMC
    ODE_SMMALA(stateData, timePoints, parameters, parmIndex, ...
                Functions, Options, SMMALAparameters);

    % Stop timer
    toc;

end