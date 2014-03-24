%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESS.m
%
% Author: Chris Burke
% Last Modified: 01-29-14
%
% Calculate ESS of MCMC sample
%
% Inputs
%   parameterHistory    Array of sampled parameter values
%   index               Index of parameter of interest
%
% Outputs
%   ess                 Effective Sample Size
%   trueVariance        Adjusted Variance of MCMC estimator
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ess, trueVariance] = ESS(parameterHistory, index)

    theta = parameterHistory(:, index);
    T = numel(theta);
    thetaMean = mean(theta);
    thetaSSE = sum((theta - thetaMean).^2);
    thetaAC = autocorr(theta, floor(T/100));
    kappa = 1 + 2*sum(thetaAC);
    
    % Efffective Sample Size
    ess = T/kappa;
    
    % Variance of MCMC estimator assuming iid sample of size ess
    trueVariance = (1/(T*ess))*thetaSSE;

end