%
% Shifted Sphere Function (S2)
% (xmin, xmax) = (-100, 100)
% error tolerance: 1e-3
%

function [ y ] = SphereS2( x )

    dim = numel(x);
    
    % x*
    o = [-25, 11, -14, 6, 50];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    z = x - o;
    
    % f(x*)
    y = -30;
    
    for i = 1 :dim
        y = y + z(i)^2;
    end
end

