%
% Shifted Rotated Griewank's Function (G3)
% (xmin, xmax) = (-600, 600)
% error tolerance: 1
%

function [ y ] = GriewankG3(x)

    dim = numel(x);
    
    % x*
    o = [11, -43, 0, 0, -540];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % linear transformation matrix, condition number 3
    T = [-0.1033   -1.0513   -0.6196    0.4665    1.2311
        -0.7389    0.8873   -0.1235   -0.6405    1.3906
        0.8542   -0.6928    1.9645    0.4103    0.1427
        -0.1204    0.4657    0.7436    1.5908    1.4934
        1.3711    0.8061   -0.5786   -1.3603    0.2375];
    
    z = (x - o)*T;
    
    % f(x*)
    y = 87;
    
    for i = 1 : dim
        y = y + (z(i)^2)/4000;
    end
    prod = 1;
    for i = 1 : dim
        prod = prod*cos(z(i)/sqrt(i));
    end
    y = y - prod + 1;
end

