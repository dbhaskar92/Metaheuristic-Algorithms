% 
% Implementation of KBPSO
%
% Author:
%
%   Dhananjay Bhaskar
%
%   A Discrete Binary Version of the Particle Swarm Algorithm, J. Kennedy,
%   C. Eberhart
%
% Last modified:
%   Friday, Oct 11, 2013
%

function [optimizer, optimum, timetot] = KBPSO(problemfile, generations, popsize)

    % params
    
    c1 = 2.0;       % cognitive parameter
    c2 = 2.0;       % social parameter
    omega = 1.0;    % inertia weight
    
    N = popsize;    % population size
    
    % initialization
    
    Q = dlmread(problemfile);    % symmetric rational matrix
    dim = size(Q,1);             % number of variables
   
    x = randi([0 1], N, dim);    % initialize position
    v = rand(N,dim);             % initialize velocity (probability)
    
    t = cputime;                        % time measurement
    gen = 0;                            % generation count
   
    pbest = x;                   % initialize best position 
    
    % initialize global best
    f_values = zeros(N,1);
    for i = 1 : N
        vectorx = x(i,:);
        sum = 0;
        for p = 1 : dim
            for q = 1 : dim
                sum = sum + Q(p,q)*vectorx(p)*vectorx(q);
            end
        end
        f_values(i,1) = sum;
    end
    
    pbest_values = f_values;
    
    [gbest_value, index] = max(f_values(:,1));
    gbest = x(index,:);
    
    while ( gen < generations )
        
            % update
            gen = gen + 1;
            
            for i = 1 : N
                
                v(i,:) = omega.*v(i,:) + c1.*rand(1,dim).*(pbest(i,:) - x(i,:)) + c2.*rand(1,dim).*(gbest - x(i,:));
               
                for j = 1 : dim
                    % update position using sigmoid limiting transformation
                    prob = 1.0/(1.0 + exp(-1.0*v(i,j)));
                    if (rand() <= prob)
                        x(i,j) = 1;
                    else
                        x(i,j) = 0;
                    end
                end
                
            end
            
            % update individual best positions
            for i = 1 : N
                vectorx = x(i,:);
                sum = 0;
                for p = 1 : dim
                    for q = 1 : dim
                        sum = sum + Q(p,q)*vectorx(p)*vectorx(q);
                    end
                end
                f_values(i,1) = sum;
            end
            
            for i = 1 : N
                if (f_values(i,1) > pbest_values(i,1))
                    
                    pbest(i,:) = x(i,:);
                    pbest_values(i,1) = f_values(i,1);
                    
                end
            end
            
            % update global best
            [candidate, pos] = max(f_values(:,1));
            if (candidate > gbest_value)
                gbest = x(pos,:);
                gbest_value = candidate;
            end
            
    end
    
    optimum = gbest_value;
    optimizer = gbest;
    timetot = cputime - t;
    
end