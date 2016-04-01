%
% Example : single variable maximization problem
%
% Author :
%   Dhananjay Bhaskar
%
% Last modified :
%   Tuesday, May 14, 2013
%

% single variable function
fun = @(x) (sin(10*x))^2/(1+x);     % global max at x = 0.1527
xmin = 0;
xmax = 1;
N = 20;                             % number of particles

% scheme for local search
scheme = 1;                         % scheme 1 or 2 (see paper for details)                    

% default values for params
cf = 0.729;
c1 = 2.05;
c2 = 2.05;

optimum = zeros(10);
optimizer = zeros(10);

% run the algorithm 10 times
for i = 1 : 10
    [optimum(i), optimizer(i)] = MPSO(N, cf, c1, c2, xmin, xmax, scheme, fun);
end

% plot the results
limits = [xmin xmax];
[X,Y] = fplot(fun, limits);

figure;
plot(X, Y, 'Color', 'r'); hold on;
plot(optimizer, optimum, '*'); hold off;
title('MPSO Optimization');
xlabel('x');
ylabel('f(x)');
