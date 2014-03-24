%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InitFHN_simple.m
%
% Author: Chris Burke
% Last Modified: 01-09-13
%
% Initial parameters for model
%
% Inputs
%
%   None
%
% Outputs
%
%   InitParm    Struct containing initial values for all parameters
%               in probabilistic model as well as fixed parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitParm = InitFHN_simple();

    % Logarithm of initial values of ODE that are assigned probability distributions
    InitParm.Y0     = [-1.0; ...      % log V0
                       1.0];          % log R0
                    
    % Logartithm of initial values of ODE that are fixed
    InitParm.Y0Fixed = [];               % No fixed initial values
    
    % Logarithm of ODE parameters that are assigned probability distributions           
    InitParm.Parm = [0.2; ...     % a
                     0.2; ...     % b
                     3.0];        % c
                    
    % Logarithm of ODE Parameters that are fixed                
    InitParm.ParmFixed = []; 
                         
    % Logarithm of standard deviations of state variable noise that are
    % assigned probability distributions
    InitParm.Sigma = [0.5; ...      % log sigmaA
                      0.5];         % log sigmaB
                    
    % Logarithm of standard deviations of state variable noise that are fixed                 
    InitParm.SigmaFixed = [];            % No fixed noise scales
                        
                    

end
