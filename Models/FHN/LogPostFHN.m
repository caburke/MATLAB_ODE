%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LogPostFHN.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% Calculates logarithm of the posterior distribution of the parameters given 
% the state data.
%
% Inputs
%
%   data    Array of states at observation times
%   times   Vector of times when observatiosn are made   
%   parm    Struct with initial parameter values
%   invTemp Inverse temperature used for tempered distributions
%
% Outputs
%
%   lpost   Logarithm of the posterior density
%   lpostd  Logarithm of the derivative of the posterior density
%   soln    Array of ODE solution values at observation times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [logPostDens, logPostDerivVal, soln, sens] = LogPostFHN(data, times, parm, invTemp, soln, sens)

% Check Positivity
if any(parm.Parm < 0) | any(parm.Sigma < 0)
    logPostDens = -1e300;
    logPostDerivVal = NaN;
    soln = NaN;
    sens = NaN;
    return;
end

%% Log Prior
logPriorDens = LogPriorFHN(parm);
logPriorDerivVal = LogPriorDerivFHN(parm);

%% Log Likelihood

% Calculate ODE Solution and Sensitivities if not provided
if nargin < 5
    [soln, sens] = SolveFHN(parm, times);
end

assert(ndims(sens) == 3, 'Sens is wrong dimension')

% Noise standard deviations used in likelihood calculation
sigma = CreateParmVector(parm, 'sigma', false);

% Check if dimensions match
[T N] = size(data);

% Calculate log likelihood
logLikeVal = 0;
for i = 1:N
    logLikeVal = logLikeVal - 0.5*T*log(2*pi*sigma(i)^2) - ...
                0.5*sum((data(:,i) - soln(:,i)).^2)/sigma(i)^2;
end

% Calculate derivative of log likelihood
[Nt, Ns, Np] = size(sens);

logLikeDeriv = zeros([Np, 1]);

% Gradient wrt model parameters
for i = 1:Np
    for j = 1:Ns
        logLikeDeriv(i) = logLikeDeriv(i) + sum((data(:, j) - soln(:, j)).*sens(:, j, i))/sigma(j)^2;
    end   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log Posterior
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
logPostDens = invTemp*(logLikeVal + logPriorDens);
logPostDerivVal = invTemp*(logLikeDeriv + logPriorDerivVal);

end

% Calculate logarithm of prior density
function logPriorDens = LogPriorFHN(parm)

% Make vectors of  parameters
p = CreateParmVector(parm, 'random', false);

logPriorDens =   LogNormalDens(p(1), -1.0, 0.2) + ...    % V0
                 LogNormalDens(p(2), 1.0, 0.2) + ...    % R0
                 LogGammaDens(p(3), 2, 1) + ...    % a
                 LogGammaDens(p(4), 2, 1) + ...    % b
                 LogGammaDens(p(5), 2, 1) + ...    % c
                 LogInvGammaDens(p(6)^2, 1, 1) + ...   % sigmaV squared
                 LogInvGammaDens(p(7)^2, 1, 1);        % sigmaR squared

end

% Calculate derivative of logarithm of prior density
function logPriorDerivVal = LogPriorDerivFHN(parm)

% Make vectors of  parameters
p = CreateParmVector(parm, 'randomODE', false);   

logPriorDerivVal = zeros([5, 1]);

logPriorDerivVal(1) = DerivLogNormalDens(p(1), -1.0, 0.2);      % log V0
logPriorDerivVal(2) = DerivLogNormalDens(p(2), 1.0, 0.2);      % log R0
logPriorDerivVal(3) = DerivLogGammaDens(p(3), 2, 1);      % log a
logPriorDerivVal(4) = DerivLogGammaDens(p(4), 2, 1);      % log b
logPriorDerivVal(5) = DerivLogGammaDens(p(5), 2, 1);      % log c

end