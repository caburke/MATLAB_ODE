%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GibbsSampleSigmaFHN.m
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

function proposedSigma = GibbsSampleSigmaFHN(currentSoln, stateData)

	[T, N] = size(currentSoln);
    
    % Hyperparameters
    a1 = 0.1;     % V
    b1 = 0.1;     % V
    a2 = 0.1;     % R
    b2 = 0.1;     % R
    
    % Calculate SSEs
    sse(1) = sum((currentSoln(:, 1) - stateData(:, 1)).^2);     % V SSE
    sse(2) = sum((currentSoln(:, 2) - stateData(:, 2)).^2);     % R SSE
    
    % Sample sigmas from conditional distributions of sigmas given other
    % parameters.  This assumes an inverse gamma prior was placed on the
    % scale parameters with hyperparameters given above.
    proposedSigma(1) = sqrt(1/gamrnd(a1 + T/2, 2/(b1 + sse(1)), 1));    % sigmaV
    proposedSigma(2) = sqrt(1/gamrnd(a2 + T/2, 2/(b2 + sse(2)), 1));    % sigmaR
    
    proposedSigma = transpose(proposedSigma);

end