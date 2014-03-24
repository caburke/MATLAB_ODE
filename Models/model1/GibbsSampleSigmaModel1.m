%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GibbsSampleSigmaModel1.m
%
% Author: Chris Burke
% Last Modified: 01-28-13
%
% Calculates the Metric Tensor used in geometric sampling algorithms
%
%   Inputs
%       currentSoln     Previous iteration solution to ODE
%       stateData       Data for model species
%
%   Outputs
%       proposedSigma   Proposed vector of noise scales
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function proposedSigma = GibbsSampleSigmaModel1(currentSoln, stateData)

	[T, N] = size(currentSoln);
    
    % Hyperparameters
    a1 = 0.1;     % A
    b1 = 0.1;     % A
    a2 = 0.1;     % B
    b2 = 0.1;     % B
    
    % Calculate SSEs
    sse(1) = sum((stateData(:, 1) - currentSoln(:, 1)).^2);     % A SSE
    sse(2) = sum((stateData(:, 2) - currentSoln(:, 2)).^2);     % B SSE
    
    % Sample sigmas from conditional distributions of sigmas given other
    % parameters.  This assumes an inverse gamma prior was placed on the
    % scale parameters with hyperparameters given above.
    proposedSigma(1) = sqrt(1/gamrnd(a1 + T/2, 1/(b1 + sse(1)/2), 1));    % sigmaA
    proposedSigma(2) = sqrt(1/gamrnd(a2 + T/2, 1/(b2 + sse(2)/2), 1));    % sigmaB
    
    proposedSigma = transpose(proposedSigma);

end