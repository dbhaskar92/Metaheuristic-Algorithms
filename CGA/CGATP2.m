%
% Test Problem 2
%
% Author :
%   Dhananjay Bhaskar
%
% Last modified :
%   Thursday, May 23, 2013
%

nvars = 2;

% function params
% maximum at x = (0.5,0.5)
xmin = 0;
xmax = 1;

% params
popsize = 30;           % population size
stringlength = 10;      % precision (binary string length)
threshold = 100;        % threshold
pc = 0.95;              % probability of crossover
pm = 0.01;              % probability of mutation
ctype = 0;              % single point crossover

% trial run
[opt, res, fe] = CGAn(@C1, nvars, xmin, xmax, pc, pm, ctype, stringlength, popsize, threshold)