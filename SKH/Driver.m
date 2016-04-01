%
%   Test Function: Ackley's Function
%   True minimum = 0;
%

dim = 20;
xmin = -32;
xmax = 32;

[minima, minimizer] = SKH(xmin, xmax, @Ackley, dim)

[minima, minimizer] = KH(xmin, xmax, @Ackley, dim)
