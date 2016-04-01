function [ y ] = Rastrigin(x)
    dim = numel(x);
    y = 10*dim;
    for i = 1 : dim
        y = y + x(i)^2 - 10*cos(2*pi*x(i));
    end
end

