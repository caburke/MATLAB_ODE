%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SecDerivLogNormalDens.m
%
% Author: Chris Burke
% Last Modified: 01-28-14
%
% Calculates second derivative of logarithm of normal density
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

function y = SecDerivLogNormalDens(x, mu, sigma);

    y = -sigma^-2;

end