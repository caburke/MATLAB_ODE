%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateParmStructFHN.m
%
% Author: Chris Burke
% Last Modified: 01-30-14
%
% Makes parameter struct out of row in MCMC parameter history
%
% Inputs
%   parmVec     Vector of random parameters in parameter history    
%
% Outputs
%
%   ParmStruct  Struct of parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ParmStruct = CreateParmStructFHN(parmVec);

    % Logarithm of initial values of ODE that are assigned probability distributions
    ParmStruct.Y0 = parmVec(1:2);
                    
    % Logartithm of initial values of ODE that are fixed
    ParmStruct.Y0Fixed = [];          % No fixed initial values
    
    % Logarithm of ODE parameters that are assigned probability distributions           
    ParmStruct.Parm = parmVec(3:5);
                    
    % Logarithm of ODE Parameters that are fixed                
    ParmStruct.ParmFixed = [];
                         
    % Logarithm of standard deviations of state variable noise that are
    % assigned probability distributions
    ParmStruct.Sigma = parmVec(6:7);
                    
    % Logarithm of standard deviations of state variable noise that are fixed                 
    ParmStruct.SigmaFixed = [];            % No fixed noise scales
                        
                    

end
