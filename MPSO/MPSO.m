% 
% Implementation of MPSO using local variant of PSO and RWDE for local
% search for unconstrained minimization problems
%
% Author:
%   Dhananjay Bhaskar
%   Memetic Particle Swarm Optimization, Y. G. Petalas, et al. 
%
% Last modified:
%   Sunday, April 17, 2016
%

function [optimum, gbest, fev, tElapsed] = MPSO(func_name,xmin,xmax,true_min,errgoal,cutoff_length,seed,c1,c2,rad,lbd,tmax)

	addpath('./function-definitions');

    % global params
    N = 60;
    dim = 5;
    fun = str2func(func_name);
    
    % local search params
    freq = 10;
    epsilon = 0.1;
    scheme = 2;
    lambda = lbd;
    
    % social params 
    kappa = 1;
    phi = c1 + c2;
    cf = 2*kappa;
    cf = cf/abs(2 - phi - sqrt(phi^2 - 4*phi));
    
    % initialization
    tStart = tic;
    
    rand('twister',seed);
    %s = RandStream('mt19937ar','Seed',seed);
    %RandStream.setGlobalStream(s);
    
    t = 0;
    fev = 0;
   
    x = xmin + (xmax - xmin).*rand(N,dim);    % initialize position
    v = xmin + (xmax - xmin).*rand(N,dim);    % initialize velocity
   
    pbest = x;                                % initialize best position 
    lbest = x;
    
    % initialize global best
    
    f_values = zeros(N,1);
    for i = 1 : N
        f_values(i,1) = feval(fun, x(i,:));
    end
    fev = fev + N;
    
    pbest_values = f_values;
    
    for i = 1 : N
        % find minimum in the neighbourhood
        [nbd_min, ~] = max(f_values(:,1));
        min_index = -1;
        for j = i-rad : i+rad
           k = j;
           if k < 1
               k = N + k;
           end
           if k > N
               k = k - N;
           end
           if f_values(k,1) <= nbd_min
               nbd_min = f_values(k,1);
               min_index = k;
           end
        end
        lbest(i,:) = x(min_index,:);
    end    
    
    [gbest_value, index] = min(f_values(:,1));
    gbest = x(index,:);                                             
    
    while (fev < cutoff_length && abs(gbest_value - true_min) > errgoal)
        
            % update
            t = t + 1;
            
            for i = 1 : N
                v(i,:) = cf.*(v(i,:) + c1.*rand(1,dim).*(pbest(i,:) - x(i,:)) + c2.*rand(1,dim).*(lbest(i,:) - x(i,:)));
                x(i,:) = x(i,:) + v(i,:);
            end
            
            % constrain x
            for i = 1 : N
                for j = 1 : dim
                	if (x(i,j) < xmin)
                		x(i,j) = xmin;
                	elseif (x(i,j) > xmax)
                		x(i,j) = xmax;
                	end 
					% random redistribution                
                    %	if ((x(i,j) < xmin) || (x(i,j) > xmax))
                    %	    x(i,j) = xmin + (xmax - xmin)*rand(1,1);
                    %	end
                end
            end
            
            % update individual best positions
            for i = 1 : N
                f_values(i,1) = feval(fun, x(i,:));
            end
            fev = fev + N;
            
            for i = 1 : N
                if (f_values(i,1) < pbest_values(i,1))
                    pbest(i,:) = x(i,:);
                    pbest_values(i,1) = f_values(i,1);
                end
            end
            
            % update local best positions
            for i = 1 : N
                % find minimum in the neighbourhood
                [nbd_min, ~] = max(pbest_values(:,1));
                min_index = -1;
                for j = i-rad : i+rad
                    k = j;
                    if k < 1
                        k = N + k;
                    end
                    if k > N
                        k = k - N;
                    end
                    if pbest_values(k,1) <= nbd_min
                        nbd_min = pbest_values(k,1);
                        min_index = k;
                    end
                end
                lbest(i,:) = pbest(min_index,:);
            end    
            
            % update global best
            [candidate, pos] = min(pbest_values(:,1));
            if (candidate < gbest_value)
                gbest = pbest(pos,:);
                gbest_value = candidate;
            end
            
            % local search (Schemes 1 and 2)
            if (mod(t,freq) == 0)
            
            	% perform local search on overall global best position
		        [y, y_value] = RWDE(fun, dim, lambda, tmax, gbest, gbest_value);
		        fev = fev + tmax;
		        if (y_value < gbest_value)
		            gbest = y;
		            gbest_value = y_value;
		        end

                if (scheme == 2)
                    
                    selection = rand(N,1) - epsilon;
                    index = find(selection < 0);
                    for i = 1 : numel(index)
                        % perform local search on pbest(index(i), :)
                        [y, y_value] = RWDE(fun, dim, lambda, tmax, pbest(index(i), :), pbest_values(index(i), 1));
                        fev = fev + tmax;
                        if (y_value < feval(fun, pbest(index(i), :)))
                            pbest(index(i), :) = y;
                            pbest_values(index(i), 1) = y_value;
                        end
                    end
                    
                end
                
            end
            
    end
    
    optimum = gbest_value;
    tElapsed = toc(tStart);
    
end


% Random walk with direct exploitation (local search)
function [y, y_value] = RWDE(fun, dim, lambda, tmax, x, x_value)

    t = 0;
    lambdainit = lambda;
    F = x_value;
    
    while (t <= tmax)
        
       t = t + 1;
       % generate unit vector
       z = rand(1, dim);
       z = z ./ norm(z,2);
       
       Fnew = feval(fun, x + lambda.*z);
       
       if (Fnew < F)
           x = x + lambda.*z;
           lambda = lambdainit;
           F = Fnew;
       end
       if (Fnew > F)
               lambda = lambda/2;
       end
       
    end
    y = x;
    y_value = F;
end
