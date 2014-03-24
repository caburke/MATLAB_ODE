%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SimulateFHN.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% Simulate data for model 1
%
% Inputs
%
%   t           Times where simulated observations occur
%   fileName    Name of file containing simulated data
%
% Outputs
%
%   None    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = SimulateFHN(t, fileName);

    % Set parameters to values in InitFHN script
    parm = InitFHN();
    sigma = CreateParmVector(parm, 'sigma', false);
    
    % Solve Model
    [y yS] = SolveFHN(parm, t);
    Soln.V = y(:,1);
    Soln.R = y(:,2);
    
    % Add noise to data
    Data.V = y(:,1) + sigma(1)*randn(size(Soln.V));
    Data.R = y(:,2) + sigma(2)*randn(size(Soln.R));
    
    % Save data to .mat file
    save(fileName, 't', 'Soln', 'Data', 'sigma', 'parm');

end