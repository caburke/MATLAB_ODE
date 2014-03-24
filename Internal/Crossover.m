%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crossover.m
%
% Author: Chris Burke
% Last Modified: 02-10-14
%
% Performs crossover step in population algorithms
%
% Inputs
%   currentParm         Struct containing parameters for all temperatures
%   currentLP           Vector of log posterior probabilities for each temp
%   invTemp
%   CreateParmStruct    Function that creates struct of parametesr
%   numRandParm         Number of random parameters
%
% Outputs
%   newParm         Struct of switched parameters
%   temp1, tmep2    Indices of switched temperatures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newParm, temp1, temp2] = Crossover(currentParm, currentLP, invTemp, CreateParmStruct, numRandParm)

    % Copy Parm Struct
    newParm = currentParm;
    
    % Variables used in function
    numTemp = numel(currentLP);
    index1 = 0;
    index2 = 0;
    temp1 = 0;
    temp2 = 0;
    
    % Create Vectors of Weights to choose temperatures
    numTemp = numel(currentLP);
    totalProb1 = 0.0;
    totalProb2 = 0.0;
    weightVec1 = zeros(numTemp, 1);
    weightVec2 = zeros(numTemp, 1);
    
    % Calculate First Weight Vector
    for n = 1:numTemp
        totalProb1 = totalProb1 + min(1e200, exp(currentLP{n}*invTemp{n}));
    end

    for n = 1:numTemp
        weightVec1(n) = min(1e200, exp(currentLP{n}*invTemp{n}))/totalProb1;
    end
    
    
    % Choose first temperature to crossover
    try
        temp1 = randsample(numTemp, 1, true, weightVec1);
    catch
        temp1 = randsample(numTemp, 1, true);        
    end
    
    % Create Second Weight Vector
    for n = 1:numTemp
        if ~(n==temp1)
            totalProb2 = totalProb2 + min(1e200, exp(currentLP{n}*invTemp{n}));
        end
    end
    
    for n = 1:numTemp
        if n == temp1
            weightVec2(n) = 0;
        else
            weightVec2(n) = min(1e200, exp(currentLP{n}*invTemp{n}))/totalProb2;
        end
    end
    
    
    % Choose second temperature to crossover
    try
        temp2 = randsample(numTemp, 1, true, weightVec2);
    catch
        temp2 = randsample(numTemp, 1, true);        
    end
    
    % Choose elements to switch
    randIndex = randsample(numRandParm, 2, false);
    index1 = randIndex(1);
    index2 = randIndex(2);
    if index1 < index2
        index = index1:index2;
    else
        index = index2:index1;
    end
    
    % Create parameter vectors for easy switching
    currentParmVec1 = CreateParmVector(currentParm{temp1}, 'random', false);
    currentParmVec2 = CreateParmVector(currentParm{temp2}, 'random', false);
    newParmVec1 = currentParmVec1;
    newParmVec2 = currentParmVec2;
    
    % Switch the chosen parameter sets
    newParmVec1(index) = currentParmVec2(index);
    newParmVec2(index) = currentParmVec1(index);
    
    % Change the vectors back into structs
    newParm{temp1} = CreateParmStruct(newParmVec1);
    newParm{temp2} = CreateParmStruct(newParmVec2);
    
end