%
% Canonical Genetic algorithm for unconstrained maximization problem
%
% Features:
%   Elitest strategy
%   Uniform and single point crossover
%
% Assumptions:
%   Population size (popsize) should be even
%   Fitness function should not be negative for any x
%
% Authors:
%   Dhananjay Bhaskar
%   Teaching Genetic Algorithms Using MATLAB, Y. J. Cao, Q. H. Wu
%
% Last modified: Thursday, May 23, 2013
%

function [optimum, optimizer, fev] = CGAn(fun, nvars, xmin, xmax, pc, pm, ctype, stringlength, popsize, threshold)

    % params
    ngenerations = 5000;        % number of generations
    
    fev = 0;                    % function evaluation count
    optcount = 0;               % optimum count
    
    % initialize population
    pop = initialise(popsize, stringlength, xmin, xmax, fun, nvars);
    fev = fev + popsize;
    [optimum, index] = max(pop(:,stringlength*nvars+nvars+1));
    
    % evolution
    gen = 1;
    
    while ((gen <= ngenerations) && (optcount <= threshold))
        
        % elitest strategy
        elite = pop(index,:);
        
        % crossover
        pop = pop(randperm(popsize),:);  % randomize mating pool
        for i = 1 : popsize/2
           parent1index = i;
           parent2index = popsize - (i-1);
           [pop(parent1index,:), pop(parent2index,:), fev] = crossover(pop(parent1index,:), pop(parent2index,:), pc, xmin, xmax, ctype, stringlength, fun, nvars, fev);
        end
        
        % mutation
        for j = 1 : popsize 
            [pop(j,:), fev] = mutation(pop(j,:), pm, xmin, xmax, stringlength, fun, nvars, fev);
        end
        
        % selection
        pop = roulette(pop, stringlength, popsize, nvars, elite);
        
        gen = gen + 1;
        
        [newopt, index] = max(pop(:,stringlength*nvars+nvars+1));
        if (newopt == optimum)
            optcount = optcount + 1; 
        else
            optcount = 0;
            optimum = newopt;
        end
        
    end
    
    optimizer = pop(index,stringlength*nvars+1:stringlength*nvars+nvars);
    
end


function [pop] = initialise(popsize, stringlength, a, b, fun, nvars)
    % randomly determine initial population
    pop = round(rand(popsize, stringlength*nvars+nvars+1));
    % calculate value of x
    for i = 1 : nvars
        bin2dvector = 2.^(size(pop(:,((i-1)*stringlength)+1:stringlength*i),2)-1:-1:0);
        pop(:,stringlength*nvars+i) = ((pop(:,((i-1)*stringlength)+1:stringlength*i)*(bin2dvector)')*(b-a))/(2.^stringlength-1)+a;
    end
    % calculate function value
    for i = 1 : popsize
        pop(i, stringlength*nvars+nvars+1) = feval(fun, pop(i,stringlength*nvars+1:stringlength*nvars+nvars));
    end
end


function [child1, child2, fev] = crossover(parent1, parent2, pc, a, b, ctype, stringlength, fun, nvars, fev)
    if (rand < pc)
        if (ctype == 0)
            % determine site for single point crossover
            cpoint = round(rand*(stringlength*nvars-1))+1;
            % perform crossover
            child1 = [parent1(:,1:cpoint) parent2(:,cpoint+1:stringlength*nvars)];
            child2 = [parent2(:,1:cpoint) parent1(:,cpoint+1:stringlength*nvars)];
        end
        if (ctype == 1)
            % perform uniform crossover
            for i = 1 : stringlength*nvars
               if (rand(1,1) < 0.5)
                   child1(:,i) = parent1(:,i);
                   child2(:,i) = parent2(:,i);
               else 
                   child1(:,i) = parent2(:,i);
                   child2(:,i) = parent1(:,i);
               end
            end
        end
        % update value of x after crossover
        for i = 1 : nvars
            bin2dvectorc1 = 2.^(size(child1(:,((i-1)*stringlength)+1:stringlength*i),2)-1:-1:0);
            bin2dvectorc2 = 2.^(size(child2(:,((i-1)*stringlength)+1:stringlength*i),2)-1:-1:0);
            child1(:,stringlength*nvars+i) = ((child1(:, ((i-1)*stringlength)+1:stringlength*i)*(bin2dvectorc1)')*(b-a))/(2.^stringlength-1)+a;
            child2(:,stringlength*nvars+i) = ((child2(:, ((i-1)*stringlength)+1:stringlength*i)*(bin2dvectorc2)')*(b-a))/(2.^stringlength-1)+a;
        end
        % update function value
        child1(:, stringlength*nvars+nvars+1) = feval(fun, child1(:, stringlength*nvars+1:stringlength*nvars+nvars));
        child2(:, stringlength*nvars+nvars+1) = feval(fun, child2(:, stringlength*nvars+1:stringlength*nvars+nvars));
        fev = fev + 2;
    else
        child1 = parent1;
        child2 = parent2;
    end
end


function [child, fev] = mutation(parent, pm, a, b, stringlength, fun, nvars, fev)
    if (rand < pm)
        % determine mutation site
        mpoint = round(rand*(stringlength*nvars-1))+1;
        child = parent;
        % flip the bit at mutation point
        child(mpoint) = abs(parent(mpoint)-1);
        % update value of x
        for i = 1 : nvars
            bin2dvectorc = 2.^(size(child(:,((i-1)*stringlength)+1:stringlength*i),2)-1:-1:0);
            child(:, stringlength*nvars+i) = ((child(:,((i-1)*stringlength)+1:stringlength*i)*(bin2dvectorc)')*(b-a))/(2.^stringlength-1)+a;
        end
        % update function value
        child(:, stringlength*nvars+nvars+1) = feval(fun, child(:, stringlength*nvars+1:stringlength*nvars+nvars));
        fev = fev + 1;
    else
        child = parent;
    end
end


function [newpop] = roulette(oldpop, stringlength, popsize, nvars, elite)

    % calculate total fitness
    totalfit = sum(oldpop(:,stringlength*nvars+nvars+1));
    % calculate relative fitness
    prob = oldpop(:,stringlength*nvars+nvars+1) / totalfit;
    % sort by cumulative fitness
    prob = cumsum(prob);
    % generate probabilites for roulette wheel selection
    rns = sort(rand(popsize,1));
    fitin = 1;    % index into old population
    newin = 1;    % index into new population
    
    % select new population, keeping population size constant
    newpop = zeros(popsize, stringlength*nvars+nvars+1);
    
    % elitest strategy
    newpop(newin,:) = elite;
    newin = newin + 1;
    
    while newin <= popsize 
        if (rns(newin) < prob(fitin))
            newpop(newin,:) = oldpop(fitin,:);
            newin = newin + 1;
        else
            fitin = fitin + 1;
        end
    end
    
end


