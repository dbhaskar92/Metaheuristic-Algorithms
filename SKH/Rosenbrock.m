function [ y ] = Rosenbrock(x)
    y = 0;
    for i = 1 : numel(x) - 1
        y = y + 100*((x(i+1) - x(i)^2)^2) + (x(i) - 1)^2;
    end
end

