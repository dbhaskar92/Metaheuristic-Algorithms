function [ y ] = Griewank(x)
    dim = numel(x);
    y = 1;
    for i = 1 : dim
        y = y + (x(i)^2)/4000;
    end
    prod = 1;
    for i = 1 : dim
        prod = prod*cos(x(i)/sqrt(i));
    end
    y = y - prod;
end

