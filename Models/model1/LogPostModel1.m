%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LogPostModel1.m
%
% Author: Chris Burke
% Last Modified: 01-28-13
%
% Calculates logarithm of the posterior distribution of the parameters given 
% the state data.
%
% Inputs
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
%   sens    Array of parameter sensitivities
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [logPostDens, logPostDerivVal, soln, sens] = LogPostModel1(data, times, parm, invTemp, soln, sens)

% Check Positivity
if any(parm.Parm < 0) || any(parm.Sigma < 0)
    logPostDens = NaN;
    logPostDerivVal = NaN;
    soln = NaN;
    sens = NaN;
    return;
end

%% Log Prior
logPriorDens = LogPriorModel1(parm);
logPriorDerivVal = LogPriorDerivModel1(parm);

%% Log Likelihood

% Calculate ODE Solution and Sensitivities if not provided
if nargin < 5
    [soln, sens, explodeFlag] = SolveModel1(parm, times);
    % Check if solution blew up
    if explodeFlag == 1;
        logPostDens = NaN;
        logPostDerivVal = NaN;
        soln = NaN;
        sens = NaN;
        return;
    end
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
    
% Only derivative wrt initial values and parameters
logLikeDeriv = zeros([Np, 1]);

% Gradient wrt model parameters
for i = 1:Np       
    for j = 1:Ns
        logLikeDeriv(i) = logLikeDeriv(i) + sum((data(:, j) - soln(:, j)).*sens(:, j, i))/sigma(j)^2;
    end    
end

%% Log Posterior

logPostDens = invTemp*(logLikeVal + logPriorDens);
logPostDerivVal = invTemp*(logLikeDeriv + logPriorDerivVal);

end

% Calculate logarithm of prior density
function logPriorDens = LogPriorModel1(parm)

% Make vector of logarithms of parameters
p = CreateParmVector(parm, 'random', false);

logPriorDens =   LogGammaDens(p(1), 2, 1) + ...    % A0
                 LogGammaDens(p(2), 2, 1) + ...    % B0
                 LogGammaDens(p(3), 2, 1) + ...    % nu
                 LogGammaDens(p(4), 2, 1) + ...    % k0
                 LogGammaDens(p(5), 2, 1) + ...    % k1
                 LogGammaDens(p(6), 2, 1) + ...    % k2
                 LogGammaDens(p(7), 2, 1) + ...    % k3
                 LogGammaDens(p(8), 2, 1) + ...    % k4
                 LogGammaDens(p(9), 2, 1) + ...    % Ka
                 LogGammaDens(p(10), 2, 1) + ...   % Kb
                 LogInvGammaDens(p(11)^2, 0.1, 0.1) + ...   % sigmaA
                 LogInvGammaDens(p(12)^2, 0.1, 0.1);        % sigmaB

end

% Calculate derivative of logarithm of prior density
function logPriorDerivVal = LogPriorDerivModel1(parm)

% Make vector of logarithms of parameters
p = CreateParmVector(parm, 'randomODE', false);

logPriorDerivVal = zeros([numel(p), 1]);

logPriorDerivVal(1) = DerivLogGammaDens(p(1), 2, 1);      % A0
logPriorDerivVal(2) = DerivLogGammaDens(p(2), 2, 1);      % B0
logPriorDerivVal(3) = DerivLogGammaDens(p(3), 2, 1);      % nu
logPriorDerivVal(4) = DerivLogGammaDens(p(4), 2, 1);      % k0
logPriorDerivVal(5) = DerivLogGammaDens(p(5), 2, 1);      % k1
logPriorDerivVal(6) = DerivLogGammaDens(p(6), 2, 1);      % k2
logPriorDerivVal(7) = DerivLogGammaDens(p(7), 2, 1);      % k3
logPriorDerivVal(8) = DerivLogGammaDens(p(8), 2, 1);      % k4
logPriorDerivVal(9) = DerivLogGammaDens(p(9), 2, 1);      % Ka
logPriorDerivVal(10) = DerivLogGammaDens(p(10), 2, 1);    % Kb

end