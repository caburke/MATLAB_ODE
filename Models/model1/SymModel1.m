%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SymModel1.m
%
% Author: Chris Burke
% Last Modified: 02-04-14
%
% Creates functions to be used in solve_model1.m
%
%   Observed States
%       A - frq mRNA
%       B - FRQ Protein
%
%   Unobserved States
%       None
%
%   Estimated Parameters
%       nu - p1
%       k0 - p2
%       k1 - p3
%       k2 - p4
%       k3 - p5
%       k4 - p6
%       Ka - p7
%       Kb - p8
%
%   Fixed Parameters
%       m - p9
%       n - p10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_states = 2;
n_parm = 10;
n_parm_est = 8;
n_sens = n_parm_est + n_states;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare Symbolic Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% States
A = sym('A');
B = sym('B');

% Parameters
nu = sym('nu');
k0 = sym('k0');
k1 = sym('k1');
k2 = sym('k2');
k3 = sym('k3');
k4 = sym('k4');
Ka = sym('Ka');
Kb = sym('Kb');
m = sym('m');   % Not estimated
n = sym('n');   % Not estimated

% Sensitivity Variables
SA = sym('SA', [1, n_sens]);
SB = sym('SB', [1, n_sens]);

% Concatenate Variables for use as vector inputs
X = [A; B];
p = [nu; k0; k1; k2; k3; k4; Ka; Kb; m; n];
p_est = [nu; k0; k1; k2; k3; k4; Ka; Kb];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Functions for Solving ODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Symbolic RHS
rhsA = k0 + nu/(1 + (B/Ka)^m) - k1*A;
rhsB = k2*A + k3*A*B^n/(Kb^n + B^n) - k4*B;
rhs = [rhsA; rhsB];

% RHS Functions
rhsA_fn = matlabFunction(rhsA, 'vars', {X, p}, 'outputs', {'dA'}, 'file', './sym_functions/model1_rhsA_fn.m');
rhsB_fn = matlabFunction(rhsB, 'vars', {X, p}, 'outputs', {'dB'}, 'file', './sym_functions/model1_rhsB_fn.m');
rhs_fn = matlabFunction(rhs, 'vars', {X, p}, 'outputs', {'ydot'}, 'file', './sym_functions/model1_rhs_fn.m');

% Jacobian
jac = jacobian(rhs, X);
jac_fn = matlabFunction(jac, 'vars', {X, p}, 'file', './sym_functions/model1_jac_fn.m');

% Parameter Derivatives
dfdp = jacobian(rhs, p_est);
dfdp_fn = matlabFunction(jacobian(rhs, p), 'vars', {X, p}, 'file', './sym_functions/model1_dfdp_fn.m');

% Symbolic First Order Sensitivities RHS
% Symbolic First Order Sensitivities RHS
SAdot = jac(1,:)*[SA; SB] + [zeros([1, n_states]), transpose(gradient(rhsA, p_est))];  % A
SBdot = jac(2,:)*[SA; SB] + [zeros([1, n_states]), transpose(gradient(rhsB, p_est))];  % B

% First Order Sensitivities RHS Functions
SAdot_fn = matlabFunction(SAdot, 'vars', {X, SA, SB, p}, 'outputs', {'SAdot'}, 'file', './sym_functions/model1_sens_rhsA_fn.m');
SBdot_fn = matlabFunction(SBdot, 'vars', {X, SA, SB, p}, 'outputs', {'SBdot'}, 'file', './sym_functions/model1_sens_rhsB_fn.m');