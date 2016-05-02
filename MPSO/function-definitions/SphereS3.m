%
% Shifted Sphere Function (S3)
% (xmin, xmax) = (-100, 100)
% error tolerance: 1e-3
%

function [ y ] = SphereS3( x )

    dim = numel(x);
    
    % x*
    o = [98, -2, 65, 34, 1];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    z = x - o;
    
    % f(x*)
    y = 111;
    
    for i = 1 :dim
        y = y + z(i)^2;
    end
end

