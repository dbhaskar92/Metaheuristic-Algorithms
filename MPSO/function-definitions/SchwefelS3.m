%
% Shifted Schwefel's Function (S3)
% (xmin, xmax) = (-100, 100)
% error tolerance: 1e-2
%

function [ y ] = SchwefelS3( x )

    dim = numel(x);
    
    % x*
    o = [0, -3, 54, 11, -7];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    z = x - o;
    
    % f(x*)
    y = 83;
    
    for i = 1 : dim
        tmp = 0;
        for j = 1 : i
            tmp = tmp + z(j);
        end
        y = y + tmp^2;
    end
end

