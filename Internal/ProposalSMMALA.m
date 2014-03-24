%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ProposalSMMALA.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% Propose move in parameter space for SMMALA algorithm
%
% Inputs
%   theta   Current parameter vector
%   gradL   Gradient of posterior distribution
%   eps     Stepsize for Euler_Maruyama integrator
%   G       Fisher Information Metric Tensor
%
% Outputs
%
%   thetaNew    Proposed parameter vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [thetaNew] = ProposalSMMALA(parm, gradL, eps, G)

    % Check dimension of G
    if not(ndims(G) == 2)
        error('Dimension of G is not 2');
    end

    % Calculate functions of FIM tensor
    try
        Ginv = inv(G);
    catch err
        fprintf('Determinant of G is %d \n', det(G));
        rethrow(err);
    end
    
    try
        RootGinv = chol(eps*Ginv);
    catch err
        fprintf('Determinant of G is %d \n', det(G));
        fprintf('Determinant of G inverse is %d \n', det(Ginv));
        fprintf('Eigenvalues of G are \n');
        eig(G)
        fprintf('Eigenvalues of G inverse are \n');
        eig(Ginv)
        rethrow(err);
    end
    
    theta = CreateParmVector(parm, 'randomODE', false);
    
    % Makes proposal by integrator one step of the SDE assuming the
    % Riemannian metric is constant near theta
    thetaNew = theta + Ginv*gradL*eps/2 + RootGinv'*randn(size(theta));

end