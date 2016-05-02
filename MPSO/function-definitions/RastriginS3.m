%
% Shifted Rotated Rastrigin's Function (S3)
% (xmin, xmax) = (-5, 5)
% error tolerance: 1
%

function [ y ] = RastriginS3(x)

    dim = numel(x);
    
    % x*
    o = [4, 1, 1, 1.5, 0];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    % Linear transformation matrix, condition number = 2
    T = [ -0.0376    0.8812   -0.6727   -0.8639   -0.4543
        -1.1217    0.1254    0.7599   -0.1782    1.0232
        0.5425   -0.7766   -0.8228   -0.2543    0.9734
        -0.0634    0.1699   -0.6882    1.2415    0.0290
        -0.4823   -0.6654   -0.5162   -0.2289   -1.0640];
    
    z = (x - o)*T;
    
    % f(x*)
    y = 98;
    
    for i = 1 : dim
        y = y + z(i)^2 - 10*cos(2*pi*z(i)) + 10;
    end
end

