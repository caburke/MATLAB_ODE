%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SimulateModel1.m
%
% Author: Chris Burke
% Last Modified: 01-15-13
%
% Simulate data for model 1
%
% Inputs
%
%   t           Times where simulated observations occur
%   sigma       Standard deviations of state data
%   fileName    Name of file containing simulated data
%
% Outputs
%
%   None    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = SimulateModel1(t, fileName);

    % Set parameters to values in InitModel1 script
    parm = InitModel1();
    sigma = CreateParmVector(parm, 'sigma', false);
    
    % Solve Model
    [y yS] = SolveModel1(parm, t);
    Soln.A = y(:,1);
    Soln.B = y(:,2);
    
    % Add noise to data
    Data.A = y(:,1) + sigma(1)*randn(size(Soln.A));
    Data.B = y(:,2) + sigma(2)*randn(size(Soln.B));
    
    % Save data to .mat file
    save(fileName, 't', 'Soln', 'Data', 'sigma', 'parm');

end