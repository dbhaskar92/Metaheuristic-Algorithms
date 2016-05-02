%
% Shifted Sphere Function (S1)
% (xmin, xmax) = (-100, 100)
% error tolerance: 1e-3
%

function [ y ] = SphereS1( x )

    dim = numel(x);
    
    % x*
    o = [-70, 3, 4, -5, 40];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    z = x - o;
    
    % f(x*)
    y = -450;
    
    for i = 1 : dim
        y = y + z(i)^2;
    end
end

