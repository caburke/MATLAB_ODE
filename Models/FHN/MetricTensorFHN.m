%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MetricTensorFHN.m
%
% Author: Chris Burke
% Last Modified: 01-28-13
%
% Calculates the Metric Tensor used in geometric sampling algorithms
%
%   Inputs
%       parm    Parameter struct
%       sens    Matrix of Parameter Sensitivities
%       invTemp
%
%   Outputs
%       metricTensor    Metric tensor for the Riemannian manifold   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculates Metric Tensor
function metricTensor = MetricTensorFHN(parm, sens, invTemp)

    if isnan(sens)
        error('Sensitivity matrix is fucked up.')
    end

    % Unpack parameters
    y0 = parm.Y0;
    p = parm.Parm;
    sigma = parm.Sigma;
    
    % Constants
    [numTimes, numStates, numParm] = size(sens);

    % Calculate Metric
    G = FisherInfoMetric(sigma, sens, numTimes, numStates, numParm);
    H = NegHessianLogPrior(parm);
    
    if not(all(size(G) == size(H)))
        G
        H
        sens
        error('Size of FIM and Hessian do not match')
    end
    metricTensor = invTemp*(G + H);
    minEig = min(eig(metricTensor));
    
    % If metric tensor is not positivie definite, make it so
    if minEig < 0
        metricTensor = metricTensor + 2*max(-minEig, 1e-5)*eye(size(metricTensor));
    end



end

%% Calculates Fisher Information Metric
function G = FisherInfoMetric(sigma, sens, numTimes, numStates, numParm)
    
    G = zeros(numParm);
    
    for i = 1:numParm
        for j = i:numParm
            for k = 1:numStates
                G(i, j) = G(i, j) + transpose(sens(:, k, i))*sens(:, k, j)/sigma(k)^2;
            end
        end
    end

    for i = 1:numParm
        for j = i:numParm
            G(j,i) = G(i,j);
        end
    end
    
end

%% Calculates Negative Hessian of Log Prior Density
function H = NegHessianLogPrior(parm)
    
    p = CreateParmVector(parm, 'randomODE', false);

    H(1,1) = - SecDerivLogNormalDens(p(1), -1.0, 0.3);
    H(2,2) = - SecDerivLogNormalDens(p(2), 1.0, 0.3);
    H(3,3) = - SecDerivLogGammaDens(p(3), 1.5, 2);
    H(4,4) = - SecDerivLogGammaDens(p(4), 1.5, 2);
    H(5,5) = - SecDerivLogGammaDens(p(5), 1.5, 2);

end