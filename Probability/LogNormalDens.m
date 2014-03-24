%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LogNormalDens.m
%
% Author: Chris Burke
% Last Modified: 01-09-13
%
% Calculates derivative of logarithm of Gamma density
%
%   Inputs
%       x       Value where function is evaluated
%       mu      Shape parameter
%       sigma   Scale parameter
%
%   Outputs
%       y   Function value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = DerivLogGammaDens(x, mu, sigma);

    y = -0.5*log(2*pi*sigma^2) - 0.5*(x - mu)^2/sigma^2;

end