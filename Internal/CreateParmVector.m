%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateParmVector.m
%
% Author: Chris Burke
% Last Modified: 01-15-13
%
% Concatenate vectors of specified type and possibly exponentiate
%
% Inputs
%
%   parm        Parameter Struct with values to concatenated
%   parmType    Type of parameters to be concatenated
%   expBool     Boolean variable specifying whether parameters will be
%               exponentiated
%
% Outputs
%
%   outVec  Vector with desired parameter type and values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outVec = CreateParmVector(parm, parmType, expBool);

    if strcmp(parmType, 'all');
        OutVec =    [parm.Y0;    parm.Y0Fixed; ...
                    parm.Parm;   parm.ParmFixed; ....
                    parm.Sigma;  parm.SigmaFixed];
                
    elseif strcmp(parmType, 'random');
        OutVec = [ parm.Y0; parm.Parm; parm.Sigma];
        
    elseif strcmp(parmType, 'ODE');
        OutVec = [ parm.Y0; parm.Y0Fixed; parm.Parm; parm.ParmFixed];
        
    elseif strcmp(parmType, 'randomODE');
        OutVec = [ parm.Y0; parm.Parm];
        
    elseif strcmp(parmType, 'fixed');
        OutVec = [ parm.Y0Fixed; parm.ParmFixed; parm.SigmaFixed];
            
    elseif strcmp(parmType, 'y0');
        OutVec = [  parm.Y0;     parm.Y0Fixed];
                 
    elseif strcmp(parmType, 'model');
        OutVec = [  parm.Parm;   parm.ParmFixed];
                 
    elseif strcmp(parmType, 'sigma');
        OutVec = [  parm.Sigma;  parm.SigmaFixed];
                 
    else;
        disp('Not a valid parameter type');
    end;
        
    if expBool == true;
        outVec = exp(OutVec);
    else;
        outVec = OutVec;
    end
end