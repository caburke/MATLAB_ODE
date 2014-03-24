%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcceptCrossover.m
%
% Author: Chris Burke
% Last Modified: 02-05-14
%
% Performs exchange step in population algorithms
%
% Inputs
%   
%
% Outputs
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [accept, flag] = AcceptCrossover(currentLP1, currentLP2, ...
                                            proposedLP1, proposedLP2)
                                     
    logAcceptProb = proposedLP1 + proposedLP2 - (currentLP1 + currentLP2);
    
    if logAcceptProb > 0
        accept = true;
    elseif log(rand()) < logAcceptProb
        accept = true;
    else
        accept = false;
    end     

end