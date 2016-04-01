%
% Genetic algorithm for unconstrained binary quadratic programming (UBQP)
%
% Features:
%   1-opt local search
%   Uniform crossover
%
% Authors:
%   Dhananjay Bhaskar
%   Genetic Algorithms for Binary Quadratic Programming, Peter Merz and
%   Bernd Freisleben
% 
% Input:
%   problemfile = 'm1bqp500.txt';   % file containing matrix Q
%   timelimit = 120;                % run time for algorithm in seconds
%
% Last modified: Monday, Aug 27, 2013
%

function [optimum, gen] = GA(problemfile, timelimit, localsearch)

    % Params
    popsize = 40;                   % population size
    
    % Initialization
    
    if (localsearch == 0)
        popsize = 100;
    end
    
    % Initialize population
    Q = dlmread(problemfile);           % symmetric rational matrix
    nvars = size(Q,1);                  % number of variables
    pop = randi([0 1], popsize, nvars); % population matrix
    
    newpop = zeros(popsize+0.5*popsize, nvars+1);   % extended population matrix for crossover
    
    gen = 0;                        % generation counter
    inittime = cputime;             % time counter
    
    % Perform local search
    if (localsearch == 1)
        for i = 1:popsize
            vectorx = pop(i,:);
            pop(i,:) = FastLocalSearch(vectorx, nvars, Q);
        end
    end
    
    % Evolution
    while ( (cputime - inittime) < timelimit )
        
        gen = gen + 1;
        
        % Crossover
        pop = pop(randperm(popsize),:);         % randomize
        
        for j = 1 : popsize
            newpop(j,1:nvars) = pop(j, :);      % copy population
        end
        
        % create new crossover vectors
        for k = 1 : 0.5*popsize
           parent1index = k;
           parent2index = popsize - (k-1);
           crossvector = CrossOver(pop, parent1index, parent2index, nvars);
           if (localsearch == 1)
                newpop(popsize+k, 1:nvars) = FastLocalSearch(crossvector, nvars, Q);
           else
                newpop(popsize+k, 1:nvars) = crossvector;
           end
        end
        
        % Selection
        
        % calculate objective function values
        for k = 1 : popsize+0.5*popsize
            vectorx = newpop(k,1:nvars);
            sum = 0;
            for p = 1 : nvars
                for q = 1 : nvars
                    sum = sum + Q(p,q)*vectorx(p)*vectorx(q);
                end
            end
            newpop(k,nvars+1) = sum;
        end
        
        % sort by objective function values : (\mu + \lambda)-ES selection
        newpop = sortrows(newpop, -1*(nvars+1));
        
        % substitute
        for k = 1 : popsize
           pop(k,:) = newpop(k,1:nvars); 
        end
        
        
        % Diversification (TODO)
        
        optimum = newpop(1,nvars+1)
        
    end
   
    % Return optimum (max) and optimizer
    optimum = newpop(1,nvars+1);

end


% 1-opt Fast Local Search
% Greedy and Local Search Heuristics for UBQP, Peter Merz and Bernd Freisleben  
function [res] = FastLocalSearch(vectorx, nvars, Q)
    
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

    while (true)
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
    

% Crossover Operator (Uniform)
function [res] = CrossOver(pop, p1index, p2index, nvars)
    parent1 = pop(p1index,:);
    parent2 = pop(p2index,:);
    
    res = zeros(1, nvars);
    for index = 1 : nvars
        if (parent1(index) == parent2(index))
            res(1,index) = parent1(index);
        else
            res(1,index) = randi([0 1]);
        end
    end
end


% Deprecated : use function above
function [res] = LocalSearch(vectorx, nvars, Q)    
    while true
        maxgain = 0;
        maxgainindex = 0;
        
        % find index for which gain is maximum
        for index = 1 : nvars 
            % calculate gain
            gain = Q(index,index)*(1 - 2*vectorx(index));
            sum = 0;
            for j = 1 : nvars
                if (j ~= index)
                   sum = sum + Q(j,index)*(vectorx(j)*(1 - 2*vectorx(index)));
                end    
            end
            gain = gain + 2*sum;
        
            if (gain > maxgain)
                maxgain = gain;
                maxgainindex = index;
            end
        end
        
        if (maxgain <= 0)
            break;
        else
            vectorx(maxgainindex) = 1 - vectorx(maxgainindex);
        end 
    end
    
    res = vectorx;    
end
