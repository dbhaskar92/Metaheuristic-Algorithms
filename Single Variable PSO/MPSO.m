% 
% Implementation of MPSO using global variant of PSO and RWDE for local
% search for single variable unconstrained maximization problems
%
% Author:
%   Dhananjay Bhaskar
%   Memetic Particle Swarm Optimization, Y. G. Petalas, et al. 
%
% Last modified:
%   Monday, May 13, 2013
%

function [optimum, gbest] = MPSO(N, cf, c1, c2, xmin, xmax, scheme, fun)

    % params
    ngenerations = 10;
    lambda = 1.0;
    tmax = 5;
    epsilon = 0.05;
    
    % initialization
    t = 0;
    x = xmin + (xmax - xmin).*rand(N,1);    % initialize position
    v = xmin + (xmax - xmin).*rand(N,1);    % initialize velocity
    pbest = x;                              % initialize best position 
    
    [~, index]= max(arrayfun(fun, x));
    gbest = x(index);                       % initialize global best      
    
    while t < ngenerations
        
            % update
            t = t + 1;
            v = cf.*(v + c1.*rand(N,1).*(pbest-x) + c2.*rand(N,1).*(gbest-x));
            x = x + v;
            
            % constrain x
            outofbounds = find(x<xmin);
            for i = 1 : numel(outofbounds)
                x(outofbounds(i)) = xmin;
            end
            outofbounds = find(x>xmax);
            for i = 1 : numel(outofbounds)
                x(outofbounds(i)) = xmax;
            end
            
            % update best positions
            newx = arrayfun(fun, x);
            diff = newx - arrayfun(fun, pbest);
            index = find(diff > 0);
            for i = 1 : numel(index)
                pbest(index(i)) = x(index(i));
            end
            
            % update global best
            [~, index]= max(arrayfun(fun, x));
            gbest = x(index);
            
            % local search (Schemes 1 and 2)
            
            if (scheme == 1)
                % local search on overall best position
                y = RWDE(fun, lambda, tmax, gbest);
                if (feval(fun, y) > feval(fun, gbest))
                    gbest = y;
                end
            end
            
            if (scheme == 2)
                selection = rand(N,1) - epsilon;
                index = find(selection<0);
                for i = 1 : numel(index)
                    % perform local search on pbest(index(i))
                    y = RWDE(fun, lambda, tmax, pbest(index(i)));
                    if (feval(fun, y) > feval(fun, pbest(index(i))))
                        pbest(index(i)) = y;
                    end
                end
            end    
            
    end
    
    optimum = feval(fun, gbest);
    
end


function [y] = RWDE(fun, lambda, tmax, x)
    t = 0;
    F = feval(fun, x);
    while (t <= tmax)
       t = t + 1;
       Fnew = feval(fun, x+lambda);
       if (Fnew < F)
           x = x + lambda;
           F = Fnew;
       end
       if (Fnew > F)
               lambda = lambda/2;
       end
    end
    y = x;
end
