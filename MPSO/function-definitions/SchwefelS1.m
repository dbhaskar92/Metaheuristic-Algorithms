%
% Shifted Schwefel's Function (S1)
% (xmin, xmax) = (-100, 100)
% error tolerance: 1e-2
%

function [ y ] = SchwefelS1( x )

    dim = numel(x);
    
    % x*
    o = [-10, 43, 41, -54, 21];
    
    if dim ~= numel(o)
        error('Error: Dimensions do not match');
    end
    
    z = x - o;
    
    % f(x*)
    y = -310;
    
    for i = 1 : dim
        tmp = 0;
        for j = 1 : i
            tmp = tmp + z(j);
        end
        y = y + tmp^2;
    end
end

