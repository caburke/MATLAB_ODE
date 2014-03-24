%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SolveModel1.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% Solves ODE Models for given parameters and returns solution
% and sensitivities wrt parameters.
%
% Inputs
%
%   parm    Struct of parameter values
%   times   Vector of the observation times
%
% Outputs
%
%   soln    Array of the ODE solution at the observation times
%   sens    Array of the sensitivities of the ODE solutions wrt the
%           parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [soln, sens, explodeFlag] = SolveModel1(parm, times)

y0 = CreateParmVector(parm, 'y0', false);    % Creates vector of initial values
p = CreateParmVector(parm, 'model', false);  % Creates vector of model parameters
% ---------------------
% CVODES initialization
% ---------------------

options = CVodeSetOptions('UserData', p, ...
                          'RelTol', 1.e-4,...
                          'AbsTol', 1.e-4, ...
                          'JacobianFn', @djacfn);

CVodeInit(@rhsfn, 'BDF', 'Newton', times(1), y0, options);

% ------------------
% FSA initialization
% ------------------

Neq = numel(y0);                                % Number of states
Ns = numel(parm.Parm) + numel(parm.Y0);   % Number of parameters to take sensitivities wrt
yS0 = [eye(Neq) zeros(Neq, Ns - Neq)];          % Initial sensitivities

% User-provided sensitivity RHS

FSAoptions = CVodeSensSetOptions('method','Staggered',...
                                 'ErrControl', true);
CVodeSensInit(Ns, @rhsSfn, yS0, FSAoptions);

% ------------------
% Integrator Output
% ------------------
nout = numel(times);
soln = zeros(nout, Neq);
sens = zeros(nout, Neq, Ns);

soln(1, :) = y0;
sens(1, :, :) = yS0;

iout = 2;
while iout <= nout,
  [status, t, y, yS] = CVode(times(iout),'Normal');
  % Exit program if solution blows up
  if any(abs(y) > 1e2)
    soln = zeros(nout, Neq);
    sens = zeros(nout, Neq, Ns);
    explodeFlag = 1;
    CVodeFree;
    return;
  end
  soln(iout, :) = y;
  sens(iout, :, :) = yS;
  if(status ==0)
    iout = iout+1;
  end
end

explodeFlag = 0;
si = CVodeGetStats;

% -----------
% Free memory
% -----------

CVodeFree;

return

% ===========================================================================
% Right-hand side function

function [yd, flag, new_data] = rhsfn(t, y, p)

yd = model1_rhs_fn(y, p);

flag = 0;
new_data = [];

return
% ===========================================================================
% Dense Jacobian function

function [J, flag, new_data] = djacfn(t, y, fy, p)

J = model1_jac_fn(y, p);

flag = 0;
new_data = [];

return

% ===========================================================================
% Sensitivity Equations RHS

function [ySd, flag, new_data] = rhsSfn(t, y, yd, yS, p)

SA = yS(1, :);
SB = yS(2, :);

ySd = zeros(size(yS));

ySd(1,:) = model1_sens_rhsA_fn(y, SA, SB, p);
ySd(2,:) = model1_sens_rhsB_fn(y, SA, SB, p);

flag = 0;
new_data = [];

return
