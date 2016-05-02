%
% Shifted Rotated Griewank's Function (G1)
% (xmin, xmax) = (-600, 600)
% error tolerance: 1
%

function [ y ] = GriewankG1(x)

    dim = numel(x);
    
    % x*
    o = [-550, 430, -110, -64, 1];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % linear transformation matrix, condition number 3
    T = [-0.0874   -1.4606   -0.8723   -0.9242   -0.3267
        -1.8564   -0.5790    1.6644   -0.6486    0.9962
        -1.3095    0.1804   -1.1988    0.5968    0.4965
        0.7168   -1.1254    0.1930    0.9313    0.2770
        0.1251    0.9453   -1.4472    0.0331    1.0578];
    
    z = (x - o)*T;
    
    % f(x*)
    y = -180;
    
    for i = 1 : dim
        y = y + (z(i)^2)/4000;
    end
    prod = 1;
    for i = 1 : dim
        prod = prod*cos(z(i)/sqrt(i));
    end
    y = y - prod + 1;
end

