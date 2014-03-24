%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ODE_Pop_SMMALA.m
%
% Author: Chris Burke
% Last Modified: 02-05-14
%
% Main Program for population simplified manifold MALA Algorithm
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

function ODE_Pop_SMMALA(stateData, timePoints, parameters, parmIndex, randParmIndex, ...
                    Functions, Options, SMMALAparameters);
    %% Initialize Algorithm
             
    % Calculate values to start MALA algorithm
    disp('Initializing SMMALA algorithm...');
    
    % Unpack run options and algorithm parameters
    burnin      = Options.burnin;
    thin        = Options.thin;
    invTemp     = Options.invTemp;
    numSamples  = Options.numSamples;
    numIter     = burnin + thin*numSamples;
    numTemp     = numel(invTemp);
    numParm     = numel(CreateParmVector(parameters{1}, 'all', false));
    numRandParm = numel(CreateParmVector(parameters{1}, 'random', false));
    numParmDeriv = numel(parmIndex.Y0) + numel(parmIndex.Parm);
    numTimes    = numel(timePoints);
    numStates   = numel(parmIndex.Y0) + numel(parmIndex.Y0Fixed);
    eps         = SMMALAparameters.StepSize;
    
    % Check inputs
    if numTemp < 3;
        error('Must include at least 3 temperatures')
    end
    
    
    % Initial Algorithm parameters
    currentParm = parameters;   % Initial copy of the parameter struct
    proposedParm = parameters;  % Copy of first temperature
    
    % Declare cell arrays for each temperature
    currentLP       = cell(numTemp, 1);
    proposedLP      = cell(numTemp, 1);
    currentLPDeriv  = cell(numTemp, 1);
    proposedLPDeriv = cell(numTemp, 1);
    currentSoln     = cell(numTemp, 1);
    proposedSoln    = cell(numTemp, 1);
    currentSens     = cell(numTemp, 1);
    proposedSens    = cell(numTemp, 1);
    currentSigma    = cell(numTemp, 1);
    currentG        = cell(numTemp, 1);
    proposedG       = cell(numTemp, 1);
    
    % Define elements of cell arrays
    for n = 1:numTemp
       currentLP{n}         = 0.0;
       proposedLP{n}        = 0.0;
       currentLPDeriv{n}    = zeros(numParmDeriv);
       proposedLPDeriv{n}   = zeros(numParmDeriv);
       currentSoln{n}       = zeros(numTimes, numStates);
       proposedSoln{n}      = zeros(numTimes, numStates);
       currentSens{n}       = zeros(numTimes, numStates, numParmDeriv);
       proposedSens{n}      = zeros(numTimes, numStates, numParmDeriv);
       currentSigma{n}      = zeros(numStates);
       currentG{n}          = zeros(numParmDeriv, numParmDeriv);
       proposedG{n}         = zeros(numParmDeriv, numParmDeriv);
    end
    
    % Initialize containers for all invTemperatures
    for n = 1:numTemp
        [currentLP{n}, currentLPDeriv{n}, currentSoln{n}, currentSens{n}] = ...
                            Functions.LogPosterior(stateData, timePoints, currentParm{n}, invTemp{n});
        currentSigma{n} = CreateParmVector(currentParm{n}, 'sigma', true);
        currentG{n} = Functions.MetricTensor(currentParm{n}, currentSens{n}, invTemp{n});
    end
    

    % Counters for loop
    acceptMutCtr    = 0;
    acceptCrossCtr  = 0;
    acceptExCtr     = 0;
    
    localAcceptMutCtr   = 0;
    localAcceptCrossCtr = 0;
    localAcceptExCtr    = 0;
    
    attemptMutCtr   = 0;
    attemptCrossCtr = 0;
    attemptExCtr    = 0;
    
    localAttemptMutCtr   = 0;
    localAttemptCrossCtr = 0;
    localAttemptExCtr    = 0;
    
    mutErrorCtr   = 0;
    crossErrorCtr = 0;
    
    sampleCtr   = 1;
    burninFlag  = true;

    % Containers for Results
    parameterHistory    = zeros([numTemp, numSamples, numRandParm]);
    logPostHistory      = zeros([numTemp, numSamples]);
    metricTensorHistory = zeros([numTemp, numSamples, numParmDeriv, numParmDeriv]);
    
    %%%%%%%%%%%%%%
    %% Main Loop
    %%%%%%%%%%%%%%
    disp('Entering Main Loop...');
    
    for i = 1:numIter
        
        if rand() < Options.MutateProb 
            %% Mutation

            % Randomly choose a temperature to mutate
            randTemp = randsample(numTemp, 1);

            % Propose noise variances using Gibbs sampling
            currentParm{randTemp}.Sigma = Functions.GibbsSampler(currentSoln{randTemp}, stateData);

            % Recalculate log posterior and log posterior derivative
            [currentLP{randTemp}, currentLPDeriv{randTemp}, currentSoln{randTemp}, currentSens{randTemp}] = ...
                Functions.LogPosterior(stateData, timePoints, currentParm{randTemp}, ...
                                    invTemp{randTemp}, currentSoln{randTemp}, currentSens{randTemp});
            currentG{randTemp} = Functions.MetricTensor(currentParm{randTemp}, currentSens{randTemp}, invTemp{randTemp});

            % Make proposal for initial values and model parameters using the SMMALA algorithm
            if Options.RandEps == true
                randEps = max(eps/10, eps + eps*randn());
            else
                randEps = eps;
            end
            
            proposedParmVector = ProposalSMMALA(currentParm{randTemp}, currentLPDeriv{randTemp}, randEps, currentG{randTemp});
            proposedParm{randTemp}.Y0 = proposedParmVector(parmIndex.Y0);
            proposedParm{randTemp}.Parm = proposedParmVector(parmIndex.Parm);
            proposedParm{randTemp}.Sigma = CreateParmVector(currentParm{randTemp}, 'sigma', true);

            % Calculate Posterior of Proposed parameters
            [proposedLP{randTemp}, proposedLPDeriv{randTemp}, proposedSoln{randTemp}, proposedSens{randTemp}] = ...
                            Functions.LogPosterior(stateData, timePoints, proposedParm{randTemp}, invTemp{randTemp}); 


            % Check if inputs for mutation are reasonable
            mutFlag = false;

            if (isfinite(proposedLP{randTemp}) && all(isfinite(proposedLPDeriv{randTemp})) && ...
                all(all(isfinite(proposedSoln{randTemp}))) && all(all(all(isfinite(proposedSens{randTemp})))))

                mutFlag = true;

            else
                mutErrorCtr = mutErrorCtr + 1;
            end
            
            if mutFlag == true
                proposedG{randTemp} = Functions.MetricTensor(proposedParm{randTemp}, proposedSens{randTemp}, invTemp{randTemp});

                % Accept SMMALA move according to MH probability
                accept = AcceptSMMALA(currentParm{randTemp}, currentLP{randTemp}, currentLPDeriv{randTemp}, currentG{randTemp}, ...
                                      proposedParm{randTemp}, proposedLP{randTemp}, proposedLPDeriv{randTemp}, proposedG{randTemp}, randEps);
                

                % Accept/Reject proposal
                if accept == true;
                    % Update containers
                    currentParm{randTemp}.Y0   = proposedParm{randTemp}.Y0;
                    currentParm{randTemp}.Parm = proposedParm{randTemp}.Parm;
                    currentLP{randTemp}        = proposedLP{randTemp};
                    currentLPDeriv{randTemp}   = proposedLPDeriv{randTemp};
                    currentG{randTemp}         = proposedG{randTemp};

                    % Increment Counters
                    acceptMutCtr        = acceptMutCtr + 1;
                    localAcceptMutCtr   = localAcceptMutCtr + 1;
                end
                
            else

                %disp('Log posterior function failed');
                
            end

            attemptMutCtr = attemptMutCtr + 1;
            localAttemptMutCtr = localAttemptMutCtr + 1;

        else
            %% Crossover
            [proposedParm, tempIndex1, tempIndex2] = Crossover(currentParm, currentLP, invTemp, ...
                                                        Functions.CreateParmStruct, numRandParm);
            invTemp1 = invTemp{tempIndex1};
            invTemp2 = invTemp{tempIndex2};
            
            [proposedLP1, proposedLPDeriv1, proposedSoln1, proposedSens1] = ...
                Functions.LogPosterior(stateData, timePoints, proposedParm{tempIndex1}, invTemp1);
            [proposedLP2, proposedLPDeriv2, proposedSoln2, proposedSens2] = ...
                Functions.LogPosterior(stateData, timePoints, proposedParm{tempIndex2}, invTemp2);

            % Check if inputs are acceptable for crossover move
            crossoverFlag = false;
            if (isfinite(proposedLP1) && isfinite(proposedLP2) && ...
                    all(isfinite(proposedLPDeriv1)) && all(isfinite(proposedLPDeriv2)) && ...
                    all(all(isfinite(proposedSoln1))) && all(all(isfinite(proposedSoln2))) && ...
                    all(all(all(isfinite(proposedSens1)))) && all(all(all(isfinite(proposedSens2)))))

                crossoverFlag = true;

            else
                crossErrorCtr = crossErrorCtr + 1;
            end

            % If both transition probabilities are reasonable, check if
            % crossover is accepted
            if crossoverFlag == true
                
                try
                    accept = AcceptCrossover(currentLP{tempIndex1}, currentLP{tempIndex2}, ...
                                            proposedLP1, proposedLP2);

                    if accept == true

                        currentParm{tempIndex1}     = proposedParm{tempIndex1};
                        currentParm{tempIndex2}     = proposedParm{tempIndex2};
                        currentLP{tempIndex1}       = proposedLP1;
                        currentLP{tempIndex2}       = proposedLP2;
                        currentLPDeriv{tempIndex1}  = proposedLPDeriv1;
                        currentLPDeriv{tempIndex2}  = proposedLPDeriv2;
                        currentSoln{tempIndex1}     = proposedSoln1;
                        currentSoln{tempIndex2}     = proposedSoln2;
                        currentSens{tempIndex1}     = proposedSens1;
                        currentSens{tempIndex2}     = proposedSens2;
                        currentG{tempIndex1}        = Functions.MetricTensor(proposedParm{tempIndex1}, ...
                                                        proposedSens1, invTemp1);
                        currentG{tempIndex2}        = Functions.MetricTensor(proposedParm{tempIndex2}, ...
                                                        proposedSens2, invTemp2);

                        acceptCrossCtr = acceptCrossCtr + 1;
                        localAcceptCrossCtr = localAcceptCrossCtr + 1;
                    end
                    
                catch err
                    fprintf('First proposed LP is %d \n', proposedLP1);
                    fprintf('Second proposed LP is %d \n', proposedLP2);
                    fprintf('Size of first sens array is %d \n', size(proposedSens1))
                    fprintf('Size of second sens array is %d \n', size(proposedSens2))
                    rethrow(err);
                end
            else
                disp('Crossover function failed');
            end
            
            attemptCrossCtr = attemptCrossCtr + 1;
            localAttemptCrossCtr = localAttemptCrossCtr + 1;
            
        end
        
        %% Exchange
        
        % Exhange Indices
        tempIndex1 = randsample(1:(numTemp-1), 1);
        tempIndex2 = tempIndex1 + 1;    
        
        % Accept Exchange
        accept = AcceptExchange(currentLP, invTemp, tempIndex1, tempIndex2);
        
        if accept == true
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Temporary data holders
            % The LP, LP derivative, and metric Tensor are multiplied by a
            % factor to account for the change in inverse temperatures that
            % occurs because of the exchange
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            invTemp1            = invTemp{tempIndex1};
            invTemp2            = invTemp{tempIndex2};
            exchangeParm1       = currentParm{tempIndex1};
            exchangeParm2       = currentParm{tempIndex2};
            exchangeLP1         = currentLP{tempIndex1}*invTemp2/invTemp1;
            exchangeLP2         = currentLP{tempIndex2}*invTemp1/invTemp2;
            exchangeLPDeriv1    = currentLPDeriv{tempIndex1}*invTemp2/invTemp1;
            exchangeLPDeriv2    = currentLPDeriv{tempIndex2}*invTemp1/invTemp2;
            exchangeSoln1       = currentSoln{tempIndex1};
            exchangeSoln2       = currentSoln{tempIndex2};
            exchangeSens1       = currentSens{tempIndex1};
            exchangeSens2       = currentSens{tempIndex2};
            exchangeG1          = currentG{tempIndex1}*invTemp2/invTemp1;
            exchangeG2          = currentG{tempIndex2}*invTemp1/invTemp2;
            
            % Now the quantities are switched
            currentParm{tempIndex1}     = exchangeParm2;
            currentParm{tempIndex2}     = exchangeParm1;
            currentLP{tempIndex1}       = exchangeLP1;
            currentLP{tempIndex2}       = exchangeLP2;
            currentLPDeriv{tempIndex1}  = exchangeLPDeriv1;
            currentLPDeriv{tempIndex2}  = exchangeLPDeriv2;
            currentSoln{tempIndex1}     = exchangeSoln1;
            currentSoln{tempIndex2}     = exchangeSoln2;
            currentSens{tempIndex1}     = exchangeSens1;
            currentSens{tempIndex2}     = exchangeSens2;
            currentG{tempIndex1}        = exchangeG1;
            currentG{tempIndex2}        = exchangeG2;
            
            acceptExCtr = acceptExCtr + 1;
            localAcceptExCtr = localAcceptExCtr + 1;
        end
        
        attemptExCtr = attemptExCtr + 1;
        localAttemptExCtr = localAttemptExCtr + 1;
        
        %% Monitor progress of MCMC
        if mod(i, 100) == 0;
            fprintf('Iteration %d of %d \n', i, numIter);
        end
        
        if i > burnin & burninFlag == true
            fprintf('###############################\n');
            fprintf('Burnin is complete\n');
            fprintf('###############################\n');
            burninFlag = false;
        end

        if (i > burnin) & (mod(i, 100) == 0);
            fprintf('\n#############################\n');
            % Acceptance Rate
            fprintf('Mutation acceptance rate is %g. \n', acceptMutCtr/attemptMutCtr);
            fprintf('Crossover acceptance rate is %g. \n', acceptCrossCtr/attemptCrossCtr);
            fprintf('Exchange acceptance rate is %g. \n', acceptExCtr/attemptExCtr);
            fprintf('Mutation local acceptance rate is %i out of %i. \n', localAcceptMutCtr, localAttemptMutCtr);
            fprintf('Crossover local acceptance rate is %i out of %i. \n', localAcceptCrossCtr, localAttemptCrossCtr);
            fprintf('Exchange local acceptance rate is %i out of %i. \n', localAcceptExCtr, localAttemptExCtr);            
            fprintf('Mutation error rate is %g. \n', mutErrorCtr/localAttemptCrossCtr);
            fprintf('Crossover error rate is %g. \n', crossErrorCtr/localAttemptCrossCtr);
            localAcceptMutCtr   = 0;
            localAcceptCrossCtr = 0;
            localAcceptExCtr    = 0;
            localAttemptMutCtr  = 0;
            localAttemptCrossCtr = 0;
            localAttemptExCtr   = 0;
            mutErrorCtr         = 0;
            crossErrorCtr       = 0;
            fprintf('\n');
        end
        
        if (i > burnin) && (mod(i, thin) == 0) && (mod(sampleCtr, 100) == 0);
            % Autocorrelation of Samples
            numLags = Options.MaxACLag;
            parmAC = zeros(numRandParm, numLags + 1);
            for j = 1:numRandParm
                parmAC(j,:) = autocorr(parameterHistory(1, 1:sampleCtr, j), numLags);
            end
            
            disp(['Maximum First Lag AC for the initial values is ', num2str(max(parmAC(randParmIndex.Y0, 2)))]);
            disp(['Minimum First Lag AC for the initial values is ', num2str(min(parmAC(randParmIndex.Y0, 2)))]);
            disp(['Maximum ', num2str(floor(numLags/2)), 'th Lag AC for the initial values is ', num2str(max(parmAC(randParmIndex.Y0, floor(numLags/2) + 1)))]);
            disp(['Minimum ', num2str(floor(numLags/2)), 'th Lag AC for the initial values is ', num2str(min(parmAC(randParmIndex.Y0, floor(numLags/2) + 1)))]);
            disp(['Maximum ', num2str(numLags), 'th Lag AC for the initial values is ', num2str(max(parmAC(randParmIndex.Y0, numLags + 1)))]);
            disp(['Minimum ', num2str(numLags), 'th Lag AC for the initial values is ', num2str(min(parmAC(randParmIndex.Y0, numLags + 1)))]);
            fprintf('\n');
            disp(['Maximum First Lag AC for model parameters is ', num2str(max(parmAC(randParmIndex.Parm, 2)))]);
            disp(['Minimum First Lag AC for model parameters is ', num2str(min(parmAC(randParmIndex.Parm, 2)))]);
            disp(['Maximum ', num2str(floor(numLags/2)), 'th Lag AC for model parameters is ', num2str(max(parmAC(randParmIndex.Parm, floor(numLags/2) + 1)))]);
            disp(['Minimum ', num2str(floor(numLags/2)), 'th Lag AC for model parameters is ', num2str(min(parmAC(randParmIndex.Parm, floor(numLags/2) + 1)))]);
            disp(['Maximum ', num2str(numLags), 'th Lag AC for model parameters is ', num2str(max(parmAC(randParmIndex.Parm, numLags + 1)))]);
            disp(['Minimum ', num2str(numLags), 'th Lag AC for model parameters is ', num2str(min(parmAC(randParmIndex.Parm, numLags + 1)))]);
            fprintf('\n');
            disp(['Maximum First Lag AC for noise scales is ', num2str(max(parmAC(randParmIndex.Sigma, 2)))]);
            disp(['Minimum First Lag AC for noise scales is ', num2str(min(parmAC(randParmIndex.Sigma, 2)))]);
            disp(['Maximum ', num2str(floor(numLags/2)), 'th Lag AC for noise scales is ', num2str(max(parmAC(randParmIndex.Sigma, floor(numLags/2) + 1)))]);
            disp(['Minimum ', num2str(floor(numLags/2)), 'th Lag AC for noise scales is ', num2str(min(parmAC(randParmIndex.Sigma, floor(numLags/2) + 1)))]);
            disp(['Maximum ', num2str(numLags), 'th Lag AC for noise scales is ', num2str(max(parmAC(randParmIndex.Sigma, numLags + 1)))]);
            disp(['Minimum ', num2str(numLags), 'th Lag AC for noise scales is ', num2str(min(parmAC(randParmIndex.Sigma, numLags + 1)))]);
            fprintf('\n');
            
            % Monitor Metric Tensor
            currentEig = eig(currentG{1});
            maxEig = max(currentEig);
            minEig = min(currentEig);
            disp(['Maximum eigenvalue of metric tensor is ', num2str(maxEig)])
            disp(['Minimum eigenvalue of metric tensor is ', num2str(minEig)])
            disp(['Ratio of max/min eigenvalues of metric tensor is ', num2str(maxEig/minEig)])
            
            %Monitor step sizes
            fprintf('Stepsize in iteration %d is %d \n', i, randEps);
        
            fprintf('#############################\n');
        end

        %% Save Iteration
        if ((i > burnin) & (mod(i, thin) == 0))
            % Loop over temperatures
            for n = 1:numTemp
                parameterHistory(n, sampleCtr, :) = transpose(CreateParmVector(currentParm{n}, 'random', false));
                logPostHistory(n, sampleCtr) = currentLP{n};
                metricTensorHistory(n, sampleCtr, :, :) = currentG{n};
            end

            % Sample was saved
            sampleCtr = sampleCtr + 1;
        end
 
    end
    
    %% Save Results
    disp('Saving Results...');
    SMMALAresults.ParameterHistory          = parameterHistory;
    SMMALAresults.MutationAcceptanceRate    = acceptMutCtr/attemptMutCtr;
    SMMALAresults.CrossoverAcceptanceRate   = acceptCrossCtr/attemptCrossCtr;
    SMMALAresults.ExchangeAcceptanceRate    = acceptExCtr/attemptExCtr;
    SMMALAresults.LogPostHistory            = logPostHistory;
    SMMALAresults.MetricTensorHistory       = metricTensorHistory;

    currentTime = fix(clock);
    fileName = ['ODE_Pop_SMMALA_' Options.resultsName '_' date '_' num2str(currentTime(4:6))];
    save(['./Results/' fileName], 'SMMALAresults');

    disp('Finished...');            

end
