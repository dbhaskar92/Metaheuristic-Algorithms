% 
% Implementation of MPSO using local variant of PSO and RWDE for local
% search for unconstrained minimization problems
%
% Author:
%   Dhananjay Bhaskar
%   Memetic Particle Swarm Optimization, Y. G. Petalas, et al. 
%
% Created:
%   Monday, May 20, 2013
%

function [opt, res, fev, seed] = MPSOln(N, lbd, eps, tmax, freq, xmin, xmax, scheme, f_handle, dim, err, fopt, seed)

    % params
    cf = 0.729;
    c1 = 2.05;
    c2 = 2.05;
    rad = 5;
    ngenerations = 10000;
    
    % initialization
    t = 0;
    fev = 0;
   
    [mat, seed] = rand_matrix(N, dim, seed);
    x = xmin + (xmax - xmin).*mat;    % initialize position
    [mat, seed] = rand_matrix(N, dim, seed);
    v = xmin + (xmax - xmin).*mat;    % initialize velocity
    
    pbest = x;                        % initialize best position 
    lbest = x;
    
    % initialize global best
    
    f_values = zeros(N,1);
    for i = 1 : N
        f_values(i,1) = f_handle(x(i,:));
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
    
    while (t < ngenerations && abs(gbest_value - fopt) > err)
        
            % update
            t = t + 1;
            
            for i = 1 : N
                [m1, seed] = rand_matrix(1, dim, seed);
                [m2, seed] = rand_matrix(1, dim, seed);
                v(i,:) = cf.*(v(i,:) + c1.*m1.*(pbest(i,:) - x(i,:)) + c2.*m2.*(lbest(i,:) - x(i,:)));
                x(i,:) = x(i,:) + v(i,:);
            end
            
            % constrain x
            for i = 1 : N
                for j = 1 : dim
                    if ((x(i,j) < xmin))
                        x(i,j) = xmin;
                    elseif ((x(i,j) > xmax))
                        x(i,j) = xmax;
                    end
                end
            end
            
            % update individual best positions
            for i = 1 : N
                f_values(i,1) = f_handle(x(i,:));
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
            
                if (scheme == 1)
                    
                    % local search on overall global best position
                    [y, y_value, seed] = RWDE(f_handle, dim, lbd, tmax, gbest, gbest_value, seed);
                    fev = fev + tmax;
                    if (y_value < gbest_value)
                        gbest = y;
                        gbest_value = y_value;
                    end
                    
                end

                if (scheme == 2)
                    
                    [random_selection, seed] = rand_matrix(N, 1, seed);
                    selection = random_selection - eps;
                    index = find(selection < 0);
                    for i = 1 : numel(index)
                        % perform local search on pbest(index(i), :)
                        [y, y_value, seed] = RWDE(f_handle, dim, lbd, tmax, pbest(index(i), :), pbest_values(index(i), 1), seed);
                        fev = fev + tmax;
                        if (y_value < f_handle(pbest(index(i), :)))
                            pbest(index(i), :) = y;
                            pbest_values(index(i), 1) = y_value;
                        end
                    end
                    
                end
                
            end
            
    end
    
    opt = gbest_value;
    res = gbest;
    
end


% Random walk with direct exploitation (local search)
function [y, y_value, seed] = RWDE(f_handle, dim, lambda, tmax, x, x_value, seed)

    t = 0;
    lambdainit = lambda;
    F = x_value;
    
    while (t <= tmax)
        
       t = t + 1;
       % generate unit vector
       [randn_matrix, seed] = generate_random_matrix(1, dim, seed);
       z = randn_matrix ./ norm(randn_matrix, 2);
       
       Fnew = f_handle(x + lambda.*z);
       
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


% Create random uniform matrix (similar to rand built-in function)
function [mat, seed] = rand_matrix(nrows, ncols, seed)
    
    mat = zeros(nrows, ncols);
    
    for row = 1 : nrows
        for col = 1 : ncols
            [value, seed] = r4_uni(seed);
            mat(row, col) = value;
        end
    end

end