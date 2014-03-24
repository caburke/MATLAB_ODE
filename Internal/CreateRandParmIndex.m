%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateRandParmIndex.m
%
% Author: Chris Burke
% Last Modified: 01-24-13
%
% Create Struct of Indices from Parameter Struct
%
% Inputs
%  parmStruct   Struct of parameter values
%
% Outputs
%
%   parmIndex   Struct with indices of different parameter categories
%               for parameters that are random (assigned priors)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parmIndex = CreateRandParmIndex(parmStruct);

    N1 = numel(parmStruct.Y0);
    N2 = numel(parmStruct.Parm);
    N3 = numel(parmStruct.Sigma);
    
    parmIndex.Y0 = 1:N1;
    parmIndex.Parm = (N1 + 1):(N1 + N2);
    parmIndex.Sigma = (N1 + N2 + 1):(N1 + N2 + N3);

end