%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run_ODE_SMMALA_FHN.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% User Interface for simplified manifold MALA Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Run_ODE_SMMALA_FHN()

    clear all;  % Clear workspace 
    tic;        % Start timer

    % Script to add paths and do other things to enable algorithm to run
    StartupODEInference();

    % Specify MCMC Options
    Options.burnin          = 0;
    Options.thin            = 1;
    Options.numSamples      = 200;
    Options.resultsName     = 'FHN';
    Options.MaxACLag        = 20;
    Options.Monitor         = 100;

    % Function Handles for Model
    Functions.LogPosterior = @LogPostFHN;
    Functions.GibbsSampler = @GibbsSampleSigmaFHN;
    Functions.MetricTensor = @MetricTensorFHN;
    Functions.CreateParmStruct = @CreateParmStructFHN;

    % Data
    load('./data/FHN_data.mat');
    stateData = [Data.V, Data.R];
    timePoints = t;

    % Initialize parameters
    weight = 0.5;
    parameters = RandInitFHN(weight);
    parmIndex = CreateRandParmIndex(parameters);

    % Specify MALA parameters   
    SMMALAparameters.StepSize     = 1e-1;    % Mean step size for Euler-Maruyama Integrator
    SMMALAparameters.ErrorBound   = 0.2;     % Upper limit on number of solver errors  

    % Run ODE_HMC
    ODE_SMMALA(stateData, timePoints, parameters, parmIndex, ...
                Functions, Options, SMMALAparameters);

    % Stop timer
    toc;

end