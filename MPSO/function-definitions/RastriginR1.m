%
% Shifted Rastrigin's Function (R1)
% (xmin, xmax) = (-5.12, 5.12)
% error tolerance: 1
%

function [ y ] = RastriginR1(x)

    dim = numel(x);
    
    % x*
    o = [-1, 2, 2, -2.5, 0];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    z = x - o;
    
    % f(x*)
    y = -330;
    
    for i = 1 : dim
        y = y + z(i)^2 - 10*cos(2*pi*z(i)) + 10;
    end
end

