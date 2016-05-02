%
% Shifted Rotated High Conditioned Elliptic Function (E1)
% (xmin, xmax) = (-100, 100)
% error tolerance: 500
%

function [ y ] = EllipticE1( x )

    dim = numel(x);
    
    % x*
    o = [83, 5, -12, -51, 4];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % random orthogonal matrix
    R = [0.2066   -0.2587   -0.3661    0.8083    0.3210
        0.0798    0.9238    0.1075    0.3518   -0.0699
        -0.2098   -0.2396    0.8560    0.3984   -0.0848
        -0.7552   -0.0201   -0.3477    0.2298   -0.5055
        -0.5801    0.1476    0.0275   -0.1070    0.7934];
    
    z = (x - o)*R;
    
    % f(x*)
    y = 145;
    
    for i = 1 : dim
        y = y + (z(i)^2)*(10^6)^((i-1)/(dim-1));
    end
end

