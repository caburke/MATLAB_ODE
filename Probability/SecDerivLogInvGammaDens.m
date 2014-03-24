%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SecDerivLogInvGammaDens.m
%
% Author: Chris Burke
% Last Modified: 01-28-14
%
% Calculates second derivative of logarithm of inverse gamma density
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

function y = SecDerivLogInvGammaDens(x, a, b);

    y = (a + 1)/x^2 - 2*b/x^3;

end