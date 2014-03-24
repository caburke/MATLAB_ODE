%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RandInitModel1.m
%
% Author: Chris Burke
% Last Modified: 02-11-14
%
% Create random starting values for each temperature
%
%   Inputs
%       weight          Weight for mixing assumed initial parameters with
%                       random initial parameters
%
%   Outputs
%      	randParmStruct  Struct containing random starting values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function randParmStruct = RandInitModel1(weight)

    p = zeros(12, 1);
    
    % Randomly sample some parameters
    p(1) = gamrnd(2, 1, 1);             % A0
    p(2) = gamrnd(2, 1, 1);             % B0
    p(3) = gamrnd(2, 1, 1);             % nu
    p(4) = gamrnd(2, 1, 1);             % k0
    p(5) = gamrnd(2, 1, 1);             % k1
    p(6) = gamrnd(2, 1, 1);             % k2
    p(7) = gamrnd(2, 1, 1);             % k3
    p(8) = gamrnd(2, 1, 1);             % k4
    p(9) = gamrnd(2, 1, 1);             % Ka
    p(10) = gamrnd(2, 1, 1);            % Kb
    p(11) = sqrt(1/gamrnd(10, 10, 1));  % sigmaA
    p(12) = sqrt(1/gamrnd(10, 10, 1));  % sigmaB
    
    randParmStruct = CreateParmStructModel1(p);
    
        % Created weighted sample
    if weight > 0;
        detParmStruct = InitModel1;
        randParmStruct.Y0 = weight*randParmStruct.Y0 + (1-weight)*detParmStruct.Y0;
        randParmStruct.Parm = weight*randParmStruct.Parm + (1-weight)*detParmStruct.Parm;
        randParmStruct.Sigma = weight*randParmStruct.Sigma + (1-weight)*detParmStruct.Sigma;
    end
    
end
