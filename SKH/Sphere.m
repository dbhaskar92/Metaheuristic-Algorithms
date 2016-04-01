function [ y ] = Sphere( x )
    y = 0;
    for i = 1 : numel(x)
        y = y + x(i)^2;
    end
    y = sqrt(y);
end

