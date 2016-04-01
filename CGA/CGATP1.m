%
% Test Problem 1 : Sphere Function
%
% Author :
%   Dhananjay Bhaskar
%
% Last modified :
%   Thursday, May 23, 2013
%

nvars = 10;

% fitness function : maximum = 1
fun = @(x) 1.0/(1.0 + ((norm(x,2))^2));
          
xmin = -100;
xmax = 100;

% params
popsize = 30;           % population size
stringlength = 20;      % precision (binary string length)
threshold = 100;        % threshold
pc = 0.80;              % probability of crossover
pm = 0.05;              % probability of mutation
ctype = 1;              % single point crossover

% trial run
[opt, res, fe] = CGAn(fun, nvars, xmin, xmax, pc, pm, ctype, stringlength, popsize, threshold)