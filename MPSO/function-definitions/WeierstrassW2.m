%
% Shifted Rotated Weierstrass Function (W2)
% (xmin, xmax) = (-0.5, 0.5)
% error tolerance: 1
%

function [ y ] = WeierstrassW2(x)

    dim = numel(x);
    a = 0.5;
    b = 3;
    kmax = 20;
    
    % x*
    o = [0.3, 0.27, 0.35, -0.4, -0.4];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % Linear transformation matrix, condition number = 5
    T = [ 0.3261   -1.5464    1.4812   -0.8294   -1.9211
        -1.8879    1.3615   -0.4180   -2.5671   -0.9174
        0.7596    0.2473   -0.2717   -0.7461    0.2046
        1.3740    2.1109    2.0480    1.2860   -0.6088
        1.5144    0.9721   -3.7200    1.3445    0.3846];
    
    z = (x - o)*T;
    
    % f(x*)
    y = -83;
    
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

