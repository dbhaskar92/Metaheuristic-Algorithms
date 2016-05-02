%
% Shifted Rotated Weierstrass Function (W1)
% (xmin, xmax) = (-0.5, 0.5)
% error tolerance: 1
%

function [ y ] = WeierstrassW1(x)

    dim = numel(x);
    a = 0.5;
    b = 3;
    kmax = 20;
    
    % x*
    o = [0.23, 0.2, -0.3, -0.1, -0.4];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % Linear transformation matrix, condition number = 5
    T = [-1.2270    0.0244    1.5381   -2.5125    0.4266
        -0.1887    2.4155   -2.6657    0.4307    1.9020
        0.7011   -2.7595   -1.1536   -1.5981    1.7597
        0.9392    2.8574   -0.0691   -0.2238   -0.8558
        -0.3098    0.4576    1.0038    1.0584    0.7063];
    
    z = (x - o)*T;
    
    % f(x*)
    y = 90;
    
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

