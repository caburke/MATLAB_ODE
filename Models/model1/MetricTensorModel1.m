%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MetricTensorModel1.m
%
% Author: Chris Burke
% Last Modified: 02-06-14
%
% Calculates the Metric Tensor used in geometric sampling algorithms
%
%   Inputs
%       parm    Parameter struct
%       sens    Matrix of Parameter Sensitivities
%       invTemp Inverse temperature used for tempered distributions
%
%   Outputs
%       G       Fisher Information Matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculates Metric Tensor
function metricTensor = MetricTensorModel1(parm, sens, invTemp)
    
    sigma = CreateParmVector(parm, 'sigma', false);

    % Constants
    [numTimes, numStates, numParm] = size(sens);
        
    % Calculate Metric Tensor
    G = FisherInfoMetric(sigma, sens, numStates, numParm);
    H = NegHessianLogPrior(parm);
    
    if not(all(size(G) == size(H)))
        G
        H
        size(sens)
        error('Size of FIM and Hessian do not match');
    end
    
    metricTensor = invTemp*(G + H);
    minEig = min(eig(metricTensor));
    
    % If metric tensor is not positivie definite, make it so
    if minEig < 0
        metricTensor = metricTensor + 100*max(-minEig, 1e-2)*eye(size(metricTensor));
    end
    
    MTInv = inv(metricTensor);
    if det(MTInv) < 1e-10
        MTinv = MTInv + 5e-1*eye(size(MTInv));
        metricTensor = inv(MTInv);
    end

end

%% Calculates Fisher Information Metric
function G = FisherInfoMetric(sigma, sens, numStates, numParm)
    
% Calculate full metric tensor
%     G = zeros(numParm);
%     
%     for i = 1:numParm
%         for j = 1:numParm
%             for k = 1:numStates
%                 G(i, j) = G(i, j) + sens(:, k, i)'*sens(:, k, j)/sigma(k)^2;
%                 if isnan(G(i,j))
%                    sens(:,k,j) 
%                 end
%             end
%         end
%     end

% Calculate diagonal tensor
    G = zeros(numParm);

    for i = 1:numParm
        for k = 1:numStates
            G(i, i) = G(i, i) + sens(:, k, i)'*sens(:, k, i)/sigma(k)^2;
            if isnan(G(i,i))
               sens(:,k,i) 
            end
        end
    end

% test tensor
%     G = eye(numParm);
    
    if isnan(G)
        sigma.^2
        size(sens(:,1,1))
        size(sens)
    end
end

%% Calculates Negative Hessian of Log Prior Density
function H = NegHessianLogPrior(parm)

    p = CreateParmVector(parm, 'randomODE', false);
    
    H(1,1) = - SecDerivLogGammaDens(p(1), 2, 1);          % A0
    H(2,2) = - SecDerivLogGammaDens(p(2), 2, 1);          % B0
    H(3,3) = - SecDerivLogGammaDens(p(3), 2, 1);          % nu
    H(4,4) = - SecDerivLogGammaDens(p(4), 2, 1);          % k0
    H(5,5) = - SecDerivLogGammaDens(p(5), 2, 1);          % k1
    H(6,6) = - SecDerivLogGammaDens(p(6), 2, 1);          % k2
    H(7,7) = - SecDerivLogGammaDens(p(7), 2, 1);          % k3
    H(8,8) = - SecDerivLogGammaDens(p(8), 2, 1);          % k4
    H(9,9) = - SecDerivLogGammaDens(p(9), 2, 1);          % Ka
    H(10,10) = - SecDerivLogGammaDens(p(10), 2, 1);       % Kb

end