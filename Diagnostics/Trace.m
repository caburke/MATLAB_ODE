%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trace.m
%
% Author: Chris Burke
% Last Modified: 01-29-14
%
% Calculate ESS of MCMC sample
%
% Inputs
%   parameterHistory    Array with parameter values
%   index               Index of variable of interest
%
% Outputs
%   none
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Trace(parameterHistory, index)

    theta = parameterHistory(:, index);
    numIter = numel(theta);
    iter = 1:numIter;
    thetaCumSum = cumsum(theta);
    size(thetaCumSum)
    size(iter)
    thetaTrace = transpose(thetaCumSum)./iter;
    
    figure;
    plot(iter, thetaTrace, '-');
    title('Trace Plot');
    
end

