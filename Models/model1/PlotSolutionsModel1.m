%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlotSolutionsModel1.m
%
% Author: Chris Burke
% Last Modified: 02-04-14
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

function [] = PlotSolutionsModel1(parameterHistory, data, origSoln)

    % Initialize Variables
    stateData = [data.A, data.B];
    trueState = [origSoln.A, origSoln.B];
    [numTimes, numStates] = size(stateData);
    [numIter, numParm] = size(parameterHistory);
    plotIter = 1:100:numIter;
    times = 1:numTimes;
    
    figure;
    for i = plotIter
       parm = CreateParmStructModel1(parameterHistory(i,:));
       [soln, sens] = SolveModel1(parm, times);
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
    title('A');
    plot(times, stateData(:, 1), '.r');
    hold on;
    subplot(2, 1, 2);
    title('B');
    plot(times, stateData(:, 2), '.b');
    
    hold off;

end