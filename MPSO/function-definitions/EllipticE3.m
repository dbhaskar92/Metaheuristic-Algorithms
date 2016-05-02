%
% Shifted Rotated High Conditioned Elliptic Function (E3)
% (xmin, xmax) = (-100, 100)
% error tolerance: 500
%

function [ y ] = EllipticE3( x )

    dim = numel(x);
    
    % x*
    o = [73, 62, 92, 26, 48];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % random orthogonal matrix
    R = [0.1935   -0.1900   -0.8718   -0.3505   -0.2088
        -0.4858    0.0707    0.1168   -0.0853   -0.8591
        -0.1017    0.4857   -0.4242    0.7567   -0.0353
        -0.3219    0.7211   -0.0803   -0.5379    0.2839
        0.7827    0.4506    0.1999   -0.0892   -0.3695];
    
    z = (x - o)*R;
    
    % f(x*)
    y = 82;
    
    for i = 1 : dim
        y = y + (z(i)^2)*(10^6)^((i-1)/(dim-1));
    end
end

