%
% Example : single variable maximization problem
%
% Author :
%   Dhananjay Bhaskar
%
% Last modified :
%   Tuesday, May 14, 2013
%

% function definition and parameters
fun = @(x) (sin(10*x))^2/(1+x);     % global max at x = 0.1527
xmin = 0;
xmax = 1;

precision = 20;                     % binary string length
ngenerations = 100;                 % number of iterations
npop = 40;                          % population size

optimum = zeros(10);
optimizer = zeros(10);

% run the algorithm 10 times
for i = 1 : 10
    [optimum(i), optimizer(i)] = CGA(fun, xmin, xmax, precision, ngenerations, npop);
end

% plot the results
limits = [xmin xmax];
[X,Y] = fplot(fun, limits);

figure;
plot(X, Y, 'Color', 'r'); hold on;
plot(optimizer, optimum, '*'); hold off;
title('GA Optimization');
xlabel('x');
ylabel('f(x)');
