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

    p = zeros(7, 1);
    
    % Randomly sample some parameters
    p(1) = -1 + 0.2*randn();            % V0
    p(2) =  1 + 0.2*randn();            % R0
    p(3) = gamrnd(2, 1, 1);             % a
    p(4) = gamrnd(2, 1, 1);             % b
    p(5) = gamrnd(2, 1, 1);             % c
    p(11) = sqrt(1/gamrnd(10, 10, 1));  % sigmaV
    p(12) = sqrt(1/gamrnd(10, 10, 1));  % sigmaR
    
    randParmStruct = CreateParmStructFHN(p);
    
        % Created weighted sample
    if weight > 0;
        detParmStruct = InitFHN();
        randParmStruct.Y0 = weight*randParmStruct.Y0 + (1-weight)*detParmStruct.Y0;
        randParmStruct.Parm = weight*randParmStruct.Parm + (1-weight)*detParmStruct.Parm;
        randParmStruct.Sigma = weight*randParmStruct.Sigma + (1-weight)*detParmStruct.Sigma;
    end
    
end
