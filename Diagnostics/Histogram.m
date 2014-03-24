%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram.m
%
% Author: Chris Burke
% Last Modified: 01-30-14
%
% Plot histogram of sample
%
% Inputs
%   parameterHistory    Array with parameter values
%   index               Index of variable of interest
%   nbin                Number of bins in histogram
%
% Outputs
%   none
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Histogram(parameterHistory, index, varname, nbin)

    % Check Inputs
    if nargin < 2
        error('Must input parameter Array and index');
    elseif nargin < 3
        varname = ['Variable ', num2str(index)];
        nbin = 25;
    elseif nargin < 4
        nbin = 25;
    end
    
    % Parameter to be plotted
    x = parameterHistory(:, index);
    
    % Histogram values
    [n, xout] = hist(x, nbin);
    
    % Kernel Density Values
    [f, xi] = ksdensity(x);
    
    figure;
    histfit(x, nbin);
    %bar(xout, n/sum(n), 'r', 'BarWidth', 1.0, 'EdgeColor', 'r');
%     hold on;
%     plot(xi, f, '-');
%     hold off;
    title(['Histogram of ', varname]);

end