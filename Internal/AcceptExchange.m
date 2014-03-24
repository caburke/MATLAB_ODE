%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AcceptExchange.m
%
% Author: Chris Burke
% Last Modified: 02-05-14
%
% Determines whether exchange step is performed
%
% Inputs
%   currentLP               Vector of current unnormalized log posterior densities
%   tempIndex1, tempIndex2  Indices to perform exchange with
%   invTemp                 Vector of inverse temperatures
%
% Outputs
%   accept      Boolean. True if exchange is accepted
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function accept = AcceptExchange(currentLP, invTemp, tempIndex1, tempIndex2)

    currentLP1 = currentLP{tempIndex1};
    currentLP2 = currentLP{tempIndex2};
    currentLP = currentLP1 + currentLP2;
    proposedLP1 = currentLP1*invTemp{tempIndex2}/invTemp{tempIndex1};
    proposedLP2 = currentLP2*invTemp{tempIndex1}/invTemp{tempIndex2};
    proposedLP = proposedLP1 + proposedLP2;

    logAcceptProb = proposedLP - currentLP;
    
    if logAcceptProb > 0
        accept = true;
    elseif log(rand()) < logAcceptProb
        accept = true;
    else
        accept = false;
    end     

end