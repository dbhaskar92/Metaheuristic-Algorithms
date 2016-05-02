%
% Shifted Rotated Rastrigin's Function (S1)
% (xmin, xmax) = (-5, 5)
% error tolerance: 1
%

function [ y ] = RastriginS1(x)

    dim = numel(x);
    
    % x*
    o = [1, 2, 3, -2.5, -2.5];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % Linear transformation matrix, condition number = 2
    T = [0.4428   -0.3414   -1.1770   -0.0940    1.3620
        -0.9983    0.1849   -0.3317    0.6957    0.5402
        0.0275    0.6405    0.7858   -0.9382    0.8316
        -0.2790    1.0951   -0.7643   -0.3059   -0.5815
        -0.3775   -0.5129    0.2124   -0.9813   -0.2991];
    
    z = (x - o)*T;
    
    % f(x*)
    y = -211;
    
    for i = 1 : dim
        y = y + z(i)^2 - 10*cos(2*pi*z(i)) + 10;
    end
end

