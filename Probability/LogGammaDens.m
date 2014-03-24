%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LogGammaDens.m
%
% Author: Chris Burke
% Last Modified: 01-09-13
%
% Calculates logarithm of Gamma density
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

function y = LogGammaDens(x, a, b);

    y = -log(gamma(a)) - a*log(b) + (a - 1)*log(x) - x/b;

end