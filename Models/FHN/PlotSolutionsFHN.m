%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlotSolutionsFHN.m
%
% Author: Chris Burke
% Last Modified: 01-30-13
%
% Plot solutions for model 1B given sampled parameters
%
% Inputs
%
%   parameterHistory    Array of sampled parameter values
%   data                Data for model species
%
% Outputs
%
%   None    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = PlotSolutionsFHN(parameterHistory, data, origSoln)

    % Initialize Variables
    stateData = [data.V, data.R];
    trueState = [origSoln.V, origSoln.R];
    [numTimes, numStates] = size(stateData);
    [numIter, numParm] = size(parameterHistory);
    plotIter = 1:100:numIter;
    times = 1:numTimes;
    
    figure;
    for i = plotIter
       parm = CreateParmStructFHN(parameterHistory(i, :));
       [soln, sens] = SolveFHN(parm, times);
       subplot(2, 1, 1)
       plot(times, soln(:, 1), '-r');
       hold on;
       subplot(2, 1, 2)
       plot(times, soln(:, 2), '-b');
       hold on;
    end
    
    subplot(2, 1, 1);
    plot(times, trueState(:, 1), '--r');
    hold on;
    subplot(2, 1, 2);
    plot(times, trueState(:, 2), '--b');
    hold on;
    
    subplot(2, 1, 1);
    title('V');
    plot(times, stateData(:, 1), '.r');
    hold on;
    subplot(2, 1, 2);
    title('R');
    plot(times, stateData(:, 2), '.b');
    
    hold off;

end