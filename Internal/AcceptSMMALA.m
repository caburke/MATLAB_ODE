%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcceptSMMALA.m
%
% Author: Chris Burke
% Last Modified: 02-11-14
%
% Determines whether SMMALA proposal is accepted
%
% Inputs
%
%   currentParm     Current parameter struct
%   currentLP       Log Posterior of current parameters
%   currentLderiv   Derivative of current log posterior
%   currentG        Metric tensor based on current parameter values
%   proposedParm    Current parameter struct
%   proposedLP      Log Posterior of current parameters
%   proposedLderiv  Derivative of current log posterior
%   proposedG       Metric tensor based on current parameter values
%   eps             Stepsize of SDE integrator
%
% Outputs
%
%   accept  Boolean variable for accept/rejecting proposed parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function accept = AcceptSMMALA(currentParm, currentLP, currentLPderiv, ...
                                        currentG, proposedParm, proposedLP, ...
                                        proposedLPderiv, proposedG, eps);
                                                             
    % Functions of Fisher Information Metric
    currentGinv = inv(currentG);
    proposedGinv = inv(proposedG); 
    if(det(currentGinv) < 1e-60)
        warning('current metric tensor inverse is nearly singular');
    end
    if(det(proposedGinv) < 1e-60)
        warning('proposed metric tensor inverse is nearly singular');
    end
                         
    % Vectors of parameters assumed random                     
    currentParmVec = CreateParmVector(currentParm, 'randomODE', false);
    proposedParmVec = CreateParmVector(proposedParm, 'randomODE', false);
                       
    % Means of transition densities
    mu1 = currentParmVec + 0.5*eps*currentGinv*currentLPderiv;
    mu2 = proposedParmVec + 0.5*eps*proposedGinv*proposedLPderiv;

    % Logarithms of transition probabilities
    try
        logNewGivenOldProb = - sum(log(diag(chol(eps*currentGinv)))) - ...
            (currentParmVec - proposedParmVec)'*(currentGinv/eps)*(currentParmVec - proposedParmVec);
        logOldGivenNewProb = - sum(log(diag(chol(eps*proposedGinv)))) - ...
            (currentParmVec - proposedParmVec)'*(proposedGinv/eps)*(currentParmVec - proposedParmVec);
    catch err
        fprintf('The determinant of the current inverse metric tensor is %d ', det(currentGinv));
        fprintf('The determinant of the proposed inverse metric tensor is %d ', det(proposedGinv));
        rethrow(err);
    end

    % Log acceptance probability
    logAccept = (proposedLP + logOldGivenNewProb) - (currentLP + logNewGivenOldProb);
    if ~isfinite(logAccept)
        warning('log acceptance probability is not finite');
    end
    
    % Determines whether proposal is accepted according to Metropolis rule
    if logAccept > 0
        accept = true;
    elseif logAccept > log(rand())
        accept = true;
    else
        accept = false;
    end
    
end