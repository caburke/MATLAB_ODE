%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% symFHN.m
%
% Author: Chris Burke
% Last Modified: 02-04-14
%
% Creates functions to be used in SolveFHN.m
%
%   Observed States
%       V - Voltage
%       R - 
%
%   Unobserved States
%       None
%
%   Estimated Parameters
%       a
%       b
%       c
%
%   Fixed Parameters
%       None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_states = 2;
n_parm = 3;
n_parm_est = 3;
n_sens = n_parm_est + n_states;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare Symbolic Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% States
V = sym('V');
R = sym('R');

% Parameters
a = sym('a');
b = sym('b');
c = sym('c');

% Sensitivity Variables
SV = sym('SV', [1, n_sens]);
SR = sym('SR', [1, n_sens]);

% Concatenate Variables for use as vector inputs
X = [V; R];
p = [a; b; c];
p_est = [a; b; c];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Functions for Solving ODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Symbolic RHS
rhsV = c*(V - (V^3)/3 + R);
rhsR = -(V - a + b*R)/c;
rhs = [rhsV; rhsR];

% RHS Functions
rhsV_fn = matlabFunction(rhsV, 'vars', {X, p}, 'outputs', {'dA'}, 'file', './sym_functions/FHN_rhsV_fn.m');
rhsR_fn = matlabFunction(rhsR, 'vars', {X, p}, 'outputs', {'dB'}, 'file', './sym_functions/FHN_rhsR_fn.m');
rhs_fn = matlabFunction(rhs, 'vars', {X, p}, 'outputs', {'ydot'}, 'file', './sym_functions/FHN_rhs_fn.m');

% Jacobian
jac = jacobian(rhs, X);
jac_fn = matlabFunction(jac, 'vars', {X, p}, 'file', './sym_functions/FHN_jac_fn.m');

% Parameter Derivatives
dfdp = jacobian(rhs, p_est);
dfdp_fn = matlabFunction(jacobian(rhs, p), 'vars', {X, p}, 'file', './sym_functions/FHN_dfdp_fn.m');

% Symbolic First Order Sensitivities RHS
SVdot = jac(1,:)*[SV; SR] + [zeros([1, n_states]), transpose(gradient(rhsV, p_est))];  % V
SRdot = jac(2,:)*[SV; SR] + [zeros([1, n_states]), transpose(gradient(rhsR, p_est))];  % R

% First Order Sensitivities RHS Functions
SVdot_fn = matlabFunction(SVdot, 'vars', {X, SV, SR, p}, 'outputs', {'Sdot'}, 'file', './sym_functions/FHN_sens_rhsV_fn.m');
SRdot_fn = matlabFunction(SRdot, 'vars', {X, SV, SR, p}, 'outputs', {'Sdot'}, 'file', './sym_functions/FHN_sens_rhsR_fn.m');
