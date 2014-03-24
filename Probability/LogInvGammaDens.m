%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LogInvGammaDens.m
%
% Author: Chris Burke
% Last Modified: 01-09-13
%
% Calculates logarithm of inverse gamma density
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

    y = a*log(b) - log(gamma(a)) - (a + 1)*log(x) - b/x;

end