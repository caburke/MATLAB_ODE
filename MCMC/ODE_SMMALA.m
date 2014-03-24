%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE_SMMALA.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% Main Program for simplified manifold MALA Algorithm
%
%   Inputs
%       stateData       Data for model species
%       timePoints      Time points where observations occur
%       parameters      Model parameters, initial values, and noise scales
%       Options         
%       MALAparameters
%       FileName        Name of file with results of MALA run
%
%   Outputs
%       None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ODE_SMMALA(stateData, timePoints, parameters, parmIndex, ...
                    Functions, Options, SMMALAparameters)
    %% Initialize Algorithm
             
    % Calculate values to start MALA algorithm
    disp('Initializing SMMALA algorithm...');
    
    % Unpack run options and algorithm parameters
    burnin = Options.burnin;
    thin = Options.thin;
    numSamples = Options.numSamples;
    numIter = burnin + thin*numSamples;
    eps = SMMALAparameters.StepSize;
    monitor = Options.Monitor;
    
    % Initial Algorithm parameters
    currentParm = parameters;   % Initial copy of the parameter struct
    proposedParm = parameters;  % another copy
    numParameters = numel(CreateParmVector(parameters, 'random', false));
    [currentLP, currentLPDeriv, currentSoln, currentSens] = ...
                        Functions.LogPosterior(stateData, timePoints, currentParm, 1);
    currentSigma = CreateParmVector(currentParm, 'sigma', false);
    currentG = Functions.MetricTensor(currentParm, currentSens, 1);
    
    % Counters for loop
    acceptCtr = 0;
    localAcceptCtr = 0;
    
    attemptCtr = 0;
    localAttemptCtr = 0;
    
    errorCtr = 0;
    sampleCtr = 1;
    burninFlag = true;

    % Containers for Outputting Results
    parameterHistory    = zeros([numSamples, numParameters]);
    logPostHistory      = zeros([numSamples, 1]);
    metricTensorHistory = cell(1, numSamples);
    
    %%%%%%%%%%%%%%
    %% Main Loop
    %%%%%%%%%%%%%%
    disp('Entering Main Loop...');
    
    for i = 1:numIter
        %% Propose noise variances using Gibbs sampling
        
        currentParm.Sigma = Functions.GibbsSampler(currentSoln, stateData);
        
        % Recalculate log posterior and log posterior derivative.
        % Provide soln and sens as inputs since they are not changed
        % by changes in the noise standard deviations.
        [currentLP, currentLPDeriv, currentSoln, currentSens] = ...
            Functions.LogPosterior(stateData, timePoints, currentParm, ...
                                1, currentSoln, currentSens);
        
        %% Make proposal for initial values and model parameters using the SMMALA algorithm
        
        % Randomly choose stepsize for integrator in proposal
        randEps = max(eps/10, eps + eps*randn());
        
        % Create proposed parameter struct
        proposedParmVector = ProposalSMMALA(currentParm, currentLPDeriv, randEps, currentG);
        proposedParm.Y0 = proposedParmVector(parmIndex.Y0);
        proposedParm.Parm = proposedParmVector(parmIndex.Parm);
        proposedParm.Sigma = CreateParmVector(currentParm, 'sigma', false);

        % Calculate Posterior of Proposed parameters
        [proposedLP, proposedLPDeriv, proposedSoln, proposedSens] = ...
                        Functions.LogPosterior(stateData, timePoints, proposedParm, 1);
        
        if proposedLP > -1e300
            proposedG = Functions.MetricTensor(proposedParm, proposedSens, 1);
            
            % Accept SMMALA move acccrding to MH probability
            accept = AcceptSMMALA(currentParm, currentLP, currentLPDeriv, currentG, ...
                                    proposedParm, proposedLP, proposedLPDeriv, proposedG, randEps);
                                
            % Flag indicates error was encountered
            if flag == 1
                errorCtr = errorCtr + 1;
            end

            % Proposal is accepted
            if accept == true;
                
                % Make proposed parameters the new current parameters
                currentParm.Y0 = proposedParm.Y0;
                currentParm.Parm = proposedParm.Parm;
                currentLP = proposedLP;
                currentLPDeriv = proposedLPDeriv;
                currentG = proposedG;
                
                % Proposal was accepted
                acceptCtr = acceptCtr + 1;
                localAcceptCtr = localAcceptCtr + 1;
            end;
            
            % Exit Program if too many errors occur
            if ((i > 100) & (errorCtr/i > SMMALAparameters.ErrorBound))
                disp('Maximum % of errors in SMMALA exceeded');
                break;
            end
        end
        
        % Proposal was attempted
        attemptCtr = attemptCtr + 1;
        localAttemptCtr = localAttemptCtr + 1;
        
                %% Monitor progress of MCMC
        if mod(i, 100) == 0;
            fprintf('Iteration %d of %d \n', i, numSamples);
        end
        
        if i > burnin & burninFlag == true
            fprintf('###############################\n');
            fprintf('Burnin is complete\n');
            fprintf('###############################\n');
            burninFlag = false;
        end

        if (i > burnin) & (mod(sampleCtr + 1, monitor) == 0) & (sampleCtr > 0);
            fprintf('\n#############################\n');
            % Acceptance Rate
            fprintf('Mutation acceptance rate is %g. \n', acceptCtr/attemptCtr);
            fprintf('Mutation local acceptance rate is %i out of %i. \n', localAcceptCtr, localAttemptCtr);
            
            % Reset local counters
            localAcceptCtr = 0;
            localAttemptCtr = 0;
            
            fprintf('\n');
            
            % Autocorrelation of Samples
            numLags = Options.MaxACLag;
            parmAC = zeros(numParameters, numLags + 1);
            for j = 1:numParameters
                parmAC(j,:) = autocorr(parameterHistory(1:sampleCtr,j), numLags);
            end
            
            % Display autocorrelation of current sampled parameters
            disp(['Maximum First Lag AC for the initial values is ', num2str(max(parmAC(parmIndex.Y0, 2)))]);
            disp(['Minimum First Lag AC for the initial values is ', num2str(min(parmAC(parmIndex.Y0, 2)))]);
            disp(['Maximum ', num2str(floor(numLags/2)), 'th Lag AC for the initial values is ', num2str(max(parmAC(parmIndex.Y0, floor(numLags/2) + 1)))]);
            disp(['Minimum ', num2str(floor(numLags/2)), 'th Lag AC for the initial values is ', num2str(min(parmAC(parmIndex.Y0, floor(numLags/2) + 1)))]);
            disp(['Maximum ', num2str(numLags), 'th Lag AC for the initial values is ', num2str(max(parmAC(parmIndex.Y0, numLags + 1)))]);
            disp(['Minimum ', num2str(numLags), 'th Lag AC for the initial values is ', num2str(min(parmAC(parmIndex.Y0, numLags + 1)))]);
            fprintf('\n');
            disp(['Maximum First Lag AC for model parameters is ', num2str(max(parmAC(parmIndex.Parm, 2)))]);
            disp(['Minimum First Lag AC for model parameters is ', num2str(min(parmAC(parmIndex.Parm, 2)))]);
            disp(['Maximum ', num2str(floor(numLags/2)), 'th Lag AC for model parameters is ', num2str(max(parmAC(parmIndex.Parm, floor(numLags/2) + 1)))]);
            disp(['Minimum ', num2str(floor(numLags/2)), 'th Lag AC for model parameters is ', num2str(min(parmAC(parmIndex.Parm, floor(numLags/2) + 1)))]);
            disp(['Maximum ', num2str(numLags), 'th Lag AC for model parameters is ', num2str(max(parmAC(parmIndex.Parm, numLags + 1)))]);
            disp(['Minimum ', num2str(numLags), 'th Lag AC for model parameters is ', num2str(min(parmAC(parmIndex.Parm, numLags + 1)))]);
            fprintf('\n');
            disp(['Maximum First Lag AC for noise scales is ', num2str(max(parmAC(parmIndex.Sigma, 2)))]);
            disp(['Minimum First Lag AC for noise scales is ', num2str(min(parmAC(parmIndex.Sigma, 2)))]);
            disp(['Maximum ', num2str(floor(numLags/2)), 'th Lag AC for noise scales is ', num2str(max(parmAC(parmIndex.Sigma, floor(numLags/2) + 1)))]);
            disp(['Minimum ', num2str(floor(numLags/2)), 'th Lag AC for noise scales is ', num2str(min(parmAC(parmIndex.Sigma, floor(numLags/2) + 1)))]);
            disp(['Maximum ', num2str(numLags), 'th Lag AC for noise scales is ', num2str(max(parmAC(parmIndex.Sigma, numLags + 1)))]);
            disp(['Minimum ', num2str(numLags), 'th Lag AC for noise scales is ', num2str(min(parmAC(parmIndex.Sigma, numLags + 1)))]);
            fprintf('\n');
            
            % Monitor Metric Tensor
            disp(['Maximum eigenvalue of metric tensor is ', num2str(max(eig(currentG)))])
            disp(['Minimum eigenvalue of metric tensor is ', num2str(min(eig(currentG)))])
            disp(['Ratio of max/min eigenvalues of metric tensor is ', num2str(max(eig(currentG))/min(eig(currentG)))])
            
            %Monitor step sizes
            disp(randEps);
            
            fprintf('#############################\n');
        end
        
        %% Save Iteration
        if ((i > burnin) & (mod(i, thin) == 0))
            parameterHistory(i, :) = transpose(CreateParmVector(currentParm, 'random', false));
            logPostHistory(i) = currentLP;
            metricTensorHistory{i} = currentG;
            
            % Sample was saved
            sampleCtr = sampleCtr + 1;
        end
 
    end
    
    %% Save Results
    disp('Saving Results...');
    SMMALAresults.parameterHistory    = parameterHistory;
    SMMALAresults.acceptanceRate      = acceptCtr/attemptCtr;
    SMMALAresults.logPostHistory      = logPostHistory;
    SMMALAresults.metricTensorHistory = metricTensorHistory;

    currentTime = fix(clock);
    fileName = ['ODE_SMMALA_' Options.resultsName '_' date '_' num2str(currentTime(4:6))];
    save(['./Results/' fileName], 'SMMALAresults');

    disp('Finished...');            

end