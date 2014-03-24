%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InitModel1.m
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

function InitParm = InitModel1()

    % Logarithm of initial values of ODE that are assigned probability distributions
    InitParm.Y0 = [  1.67; ...      % A0
                        3.85];      % B0
                    
    % Logartithm of initial values of ODE that are fixed
    InitParm.Y0Fixed = [];               % No fixed initial values
    
    % Logarithm of ODE parameters that are assigned probability distributions           
    InitParm.Parm = [1.177; ...     % nu
                        0.01; ...      % k0
                        0.08; ...      % k1
                        0.0482; ...    % k2
                        1.605; ...     % k3
                        0.535; ...     % k4
                        1.1; ...       % Ka
                        3.0];          % Kb
                    
    % Logarithm of ODE Parameters that are fixed                
    InitParm.ParmFixed = [3.0; ...  % m
                             2.0];  % n
                         
    % Logarithm of standard deviations of state variable noise that are
    % assigned probability distributions
    InitParm.Sigma = [0.2; ...      % sigmaA
                      0.2];      % sigmaB
                    
    % Logarithm of standard deviations of state variable noise that are fixed                 
    InitParm.SigmaFixed = [];       % No fixed noise scales
                        
                    

end
