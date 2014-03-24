%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exchange.m
%
% Author: Chris Burke
% Last Modified: 02-05-14
%
% Performs exchange step in population algorithms
%
% Inputs
%   currentParm     Struct of current parameter values
%   numTemp         Number of temperatures used
%
% Outputs
%   newParm         Struct of new parameter values
%   temp1, temp2    Indices of switched temperatures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newParm = Exchange(currentParm, numTemp)
    
    % Exhange Indices
    temp1 = randsample(numTemp);
    if temp1 == 1
        temp2 = 2
    elseif temp1 == numTemp
        temp2 = numTemp - 1;
    else
        temp2 = temp1 + randsample([-1, 1], 1);
    end
    index = min([temp1, temp2]):max([temp1, temp2]);
    
    % Exchange the parameter sets
    newParm{temp1} = currentParm{temp2};
    newParm{temp2} = currentParm{temp1};

end