%
% Shifted Rotated Rastrigin's Function (S2)
% (xmin, xmax) = (-5, 5)
% error tolerance: 1
%

function [ y ] = RastriginS2(x)

    dim = numel(x);
    
    % x*
    o = [3, -4, 3, -2.5, 0];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % Linear transformation matrix, condition number = 2
    T = [0.1808   -0.0847   -0.0553    0.7277   -1.0206
        0.0385   -0.7343   -0.2298   -0.8880   -1.2422
        0.1196    1.1420   -0.3345   -0.5567   -0.7156
        -0.7181   -0.0196    1.5177   -0.2282   -0.1149
        1.0690   -0.2992    0.9422    0.2736    0.0631];
    
    z = (x - o)*T;
    
    % f(x*)
    y = 41;
    
    for i = 1 : dim
        y = y + z(i)^2 - 10*cos(2*pi*z(i)) + 10;
    end
end

