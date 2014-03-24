%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateParmStructModel1.m
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

function ParmStruct = CreateParmStructModel1(parmVec);

    % Logarithm of initial values of ODE that are assigned probability distributions
    ParmStruct.Y0 = parmVec(1:2);
                    
    % Logartithm of initial values of ODE that are fixed
    ParmStruct.Y0Fixed = [];          % No fixed initial values
    
    % Logarithm of ODE parameters that are assigned probability distributions           
    ParmStruct.Parm = parmVec(3:10);
                    
    % Logarithm of ODE Parameters that are fixed                
    ParmStruct.ParmFixed = [3.0; 2.0];
                         
    % Logarithm of standard deviations of state variable noise that are
    % assigned probability distributions
    ParmStruct.Sigma = parmVec(11:12);
                    
    % Logarithm of standard deviations of state variable noise that are fixed                 
    ParmStruct.SigmaFixed = [];            % No fixed noise scales
                        
                    

end
