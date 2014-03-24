%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DerivLogGammaDens.m
%
% Author: Chris Burke
% Last Modified: 01-09-13
%
% Calculates derivative of logarithm of Gamma density
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

function y = DerivLogGammaDens(x, a, b);

    y = (a - 1)/x - (1/b);

end