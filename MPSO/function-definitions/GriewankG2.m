%
% Shifted Rotated Griewank's Function (G2)
% (xmin, xmax) = (-600, 600)
% error tolerance: 1
%

function [ y ] = GriewankG2(x)

    dim = numel(x);
    
    % x*
    o = [150, 70, 10, 3, -54];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % linear transformation matrix, condition number 3
    T = [-0.4967    0.2830   -1.3181    0.3170   -0.8621
        0.9844    0.0992   -0.8604    0.0229   -1.7113
        0.2426    2.0618   -0.3053    0.9080    0.7775
        0.1395   -0.6754    0.3958    1.5883   -1.3978
        0.1684   -1.4719   -0.5880    1.1218    0.4518];
    
    z = (x - o)*T;
    
    % f(x*)
    y = 4;
    
    for i = 1 : dim
        y = y + (z(i)^2)/4000;
    end
    prod = 1;
    for i = 1 : dim
        prod = prod*cos(z(i)/sqrt(i));
    end
    y = y - prod + 1;
end

