%
% Shifted Rotated Weierstrass Function (W3)
% (xmin, xmax) = (-0.5, 0.5)
% error tolerance: 1
%

function [ y ] = WeierstrassW3(x)

    dim = numel(x);
    a = 0.5;
    b = 3;
    kmax = 20;
    
    % x*
    o = [0.115, -0.12, -0.45, 0.2, 0.32];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % Linear transformation matrix, condition number = 5
    T = [1.4217   -1.0263   -0.8829    0.5062   -2.5511
        0.6480    0.3144   -0.3768   -0.6523    0.3407
        -1.2645    0.9484    2.7669   -1.9606   -1.4636
        -3.0370   -0.4535   -1.5683   -0.8020   -0.0603
        2.8447   -1.9885    1.3280   -0.5275    0.6449];
    
    z = (x - o)*T;
    
    % f(x*)
    y = 118;
    
    for i = 1 : dim
        sum1 = 0;
        for k = 0 : kmax
            sum1 = sum1 + (a^k)*cos(2*pi*(b^k)*(z(i)+0.5));
        end
        y = y + sum1;
    end
    sum2 = 0;
    for k = 0 : kmax
        sum2 = sum2 + (a^k)*cos(2*pi*(b^k)*0.5);
    end
    y = y - dim*sum2;
end

