%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DerivLogInvGammaDens.m
%
% Author: Chris Burke
% Last Modified: 01-09-13
%
% Calculates derivative of logarithm of inverse gamma density
%
%   Inputs
%       x   Value where function is evaluated
%       a   Shape parameter
%       b   Scale parameter
%
%   Outputs
%       y   Function value
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = DerivLogInvGammaDens(x, a, b);

    y = -(a + 1)/x + b/x^2;

end