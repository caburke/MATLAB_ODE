%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateParmIndex.m
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

function parmIndex = CreateParmIndex(parmStruct);

    N1 = numel(parmStruct.Y0);
    N2 = numel(parmStruct.Y0Fixed);
    N3 = numel(parmStruct.Parm);
    N4 = numel(parmStruct.ParmFixed);
    N5 = numel(parmStruct.Sigma);
    N6 = numel(parmStruct.SigmaFixed);
    
    parmIndex.Y0 = 1:N1;
    parmIndex.Y0Fixed = (N1 + 1):(N1 + N2);
    parmIndex.Parm = (N1 + N2 + 1):(N1 + N2 + N3);
    parmIndex.ParmFixed = (N1 + N2 + N3 + 1):(N1 + N2 + N3 + N4);
    parmIndex.Sigma = (N1 + N2 + N3 + N4 + 1):(N1 + N2 + N3 + N4 + N5);
    parmIndex.SigmaFixed = (N1 + N2 + + N3 + N4 + N5 + 1):(N1 + N2 + N3 + N4 + N5 + N6);

end