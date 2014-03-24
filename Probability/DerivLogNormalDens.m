%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DerivLogNormalDens.m
%
% Author: Chris Burke
% Last Modified: 01-09-13
%
% Calculates derivative of logarithm of normal density
%
%   Inputs
%       x       Value where function is evaluated
%       mu      Shape parameter
%       sigma   Scale parameter
%
%   Outputs
%       y       Function value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = DerivLogNormalDens(x, mu, sigma);

    y = -(x - mu)/sigma^2;

end