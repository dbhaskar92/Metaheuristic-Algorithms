%
% Shifted Rotated High Conditioned Elliptic Function (E2)
% (xmin, xmax) = (-100, 100)
% error tolerance: 500
%

function [ y ] = EllipticE2( x )

    dim = numel(x);
    
    % x*
    o = [-57, -64, 28, 2, -8];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % random orthogonal matrix
    R = [-0.3788   -0.4286   -0.2511   -0.4270   -0.6538
        0.1229   -0.4046   -0.7359    0.5087    0.1444
        0.3212   -0.0329    0.3382    0.5759   -0.6706
        0.0272   -0.8017    0.5090    0.0034    0.3120
        0.8588   -0.0934   -0.1481   -0.4767   -0.0681];
    
    z = (x - o)*R;
    
    % f(x*)
    y = -51;
    
    for i = 1 : dim
        y = y + (z(i)^2)*(10^6)^((i-1)/(dim-1));
    end
end

