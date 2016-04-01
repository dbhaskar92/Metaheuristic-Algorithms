% 
% Implementation of MPSO using global variant of PSO and local search
% for unconstrained binary quadratic programming (UBQP)
%
% Author:
%
%   Dhananjay Bhaskar
%
%   Particle Swarm Optimization for Integer Programming, Laskari,
%   Parsopoulos, Vrahatis
%
%   A Discrete Binary Version of the Particle Swarm Algorithm, J. Kennedy,
%   C. Eberhart
%
% Last modified:
%   Wednesday, Aug 28, 2013
%

function [optimizer, optimum, gen] = MPSO_opt(problemfile, timelimit)

    % params
    cf = 0.729;     % constriction factor
    c1 = 2.0;       % cognitive parameter
    c2 = 2.0;       % social parameter
    omega = 1.0;    % inertia weight
    N = 40;         % population size
    Vmax = 6.0;     % maximum allowed velocity
    
    scheme = 1;     % 0: disable local search, 1: enable 
    epsilon = 0.2;  % local search probability
    
    % initialization
    
    Q = dlmread(problemfile);    % symmetric rational matrix
    dim = size(Q,1);             % number of variables
   
    x = randi([0 1], N, dim);    % initialize position
    v = rand(N,dim);             % initialize velocity (probability)
    
    t = cputime;                        % time measurement
    gen = 0;                            % generation count
    slope = (0.1 - 1.0) / (timelimit);  % slope for omega
    
    % perform local search
    if (scheme == 1)
        gainMat = zeros(N,dim);
        % calculate gains for all i in {1,...,n} ) in O(n^2)
        for index = 1 : dim 
            % calculate gain
            gainMat(:,index) = Q(index,index).*(1 - 2.*x(:,index));
            sum = zeros(N,1);
            for j = 1 : dim
                if (j ~= index)
                    sum(:,1) = sum(:,1) + Q(j,index).*(x(:,j).*(1 - 2.*x(:,index)));
                end    
            end
            gainMat(:,index) = gainMat(:,index) + 2.*sum(:,1);
        end
            
        for i = 1 : N
            [x(i,:), gainMat(i,:)] = FastLocalSearch(x(i,:), dim, Q, gainMat(i,:));
        end
    end
   
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
    
    while ( (cputime - t) < timelimit )
        
            % update
            gen = gen + 1;
            
            for i = 1 : N
                v(i,:) = omega.*v(i,:) + c1.*rand(1,dim).*(pbest(i,:) - x(i,:)) + c2.*rand(1,dim).*(gbest - x(i,:));
                %v(i,:) = cf.*v(i,:);
               
                for j = 1 : dim
                    % constrain velocity
                    if ( abs(v(i,j)) > Vmax)
                        if (v(i,j) > 0)
                            v(i,j) = Vmax;
                        else
                            v(i,j) = -1*Vmax;
                        end
                    end
                    % update position using sigmoid limiting transformation
                    prob = 1.0/(1.0 + exp(-1.0*v(i,j)));
                    if (rand() < prob)
                        x(i,j) = 1;
                    else
                        x(i,j) = 0;
                    end
                end
            end
            
            
            % linearly decrease omega
            % omega = slope*(cputime-t) + 1.0;
            
            
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
            
            if (scheme == 1)
        
                % perform local search on pbest
                selection = rand(N,1) - epsilon;
                ind = find(selection < 0);
                for i = 1 : numel(ind)
                    
                    % update gain
                    for index = 1 : dim 
                        % calculate gain
                        gainMat(ind(i),index) = Q(index,index)*(1 - 2*pbest(ind(i),index));
                        sum = 0;
                        for j = 1 : dim
                            if (j ~= index)
                                sum = sum + Q(j,index)*(pbest(ind(i),j)*(1 - 2*pbest(ind(i),index)));
                            end    
                        end
                        gainMat(ind(i),index) = gainMat(ind(i),index) + 2*sum;
                    end
                    
                    % local search
                    [y, gainMat(ind(i),:)] = FastLocalSearch(pbest(ind(i),:), dim, Q, gainMat(ind(i),:));
                    
                    % update pbest and gbest
                    pbest(ind(i),:) = y;
                    sum = 0;
                    for p = 1 : dim
                        for q = 1 : dim
                            sum = sum + Q(p,q)*y(p)*y(q);
                        end
                    end
                    y_value = sum;
                    pbest_values(ind(i), 1) = y_value;
                    if (y_value > gbest_value)
                        gbest = y;
                        gbest_value = y_value;
                    end
                end
            end
            
            optimum = gbest_value
            
    end
    
    optimum = gbest_value;
    optimizer = gbest;
    
end


% Local Search
function [res, gain] = FastLocalSearch(vectorx, nvars, Q, gain)
    
    % perform local search in O(n)
        [maxgain, maxgainindex] = max(gain);
        if (maxgain > 0)
            vectorx(maxgainindex) = 1 - vectorx(maxgainindex);
            % update gains
            for index = 1 : nvars
                if (index == maxgainindex)
                    gain(index) = -1*gain(index);
                else
                    if (Q(index,maxgainindex) ~= 0)
                        product = (1 - 2*vectorx(index))*(2*vectorx(maxgainindex) - 1);
                        product = product*2*Q(index,maxgainindex);
                        gain(index) = gain(index) + product;
                    end    
                end 
            end
        end
    res = vectorx;
end


% 1-opt Fast Local Search
% Greedy and Local Search Heuristics for UBQP, Peter Merz and Bernd Freisleben  
function [res] = LocalSearchAlgorithm(vectorx, nvars, Q)
    
    gain = zeros(1,nvars);
    
    % calculate gains for all i in {1,...,n} ) in O(n^2)
    for index = 1 : nvars 
        % calculate gain
        gain(index) = Q(index,index)*(1 - 2*vectorx(index));
        sum = 0;
        for j = 1 : nvars
            if (j ~= index)
                sum = sum + Q(j,index)*(vectorx(j)*(1 - 2*vectorx(index)));
            end    
        end
        gain(index) = gain(index) + 2*sum;
    end
    
    % perform local search in O(n)
    while true
        [maxgain, maxgainindex] = max(gain);
        if (maxgain <= 0)
            break;
        else
            vectorx(maxgainindex) = 1 - vectorx(maxgainindex);
            % update gains
            for index = 1 : nvars
                if (index == maxgainindex)
                    gain(index) = -1*gain(index);
                else
                    if (Q(index,maxgainindex) ~= 0)
                        product = (1 - 2*vectorx(index))*(2*vectorx(maxgainindex) - 1);
                        product = product*2*Q(index,maxgainindex);
                        gain(index) = gain(index) + product;
                    end    
                end 
            end
        end
    end
    res = vectorx;
end