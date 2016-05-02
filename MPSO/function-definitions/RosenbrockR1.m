%
% Shifted Rosenbrocks's Function (R1)
% (xmin, xmax) = (-30, 30)
% error tolerance: 1
%

function [ y ] = RosenbrockR1(x)
    
    dim = numel(x);
    
    % x*
    o = [0, -3, 0, -25, 10];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    z = x - o + 1;
    
    % f(x*)
    y = 390;
    
    for i = 1 : dim - 1
        y = y + 100*((z(i)^2 - z(i+1))^2) + (z(i) - 1)^2;
    end
end

