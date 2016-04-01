%
% Genetic algorithm for single variable unconstrained maximization problem
%
% Assumptions:
%   fun(x) > 0 for all xmin<=x<=xmax
%   population size (popsize) should be even
%
% Authors:
%   Dhananjay Bhaskar
%   Teaching Genetic Algorithms Using MATLAB, Y. J. Cao, Q. H. Wu
%
% Last modified: Monday, May 13, 2013
%

function [optimum, optimizer] = CGA(fun, xmin, xmax, stringlength, ngenerations, popsize)

    % params
    pc = 0.95;      % probability of crossover
    pm = 0.05;      % probability of mutation
    
    % initialize population
    pop = initialise(popsize, stringlength, xmin, xmax, fun);
    
    % evolution
    gen = 1;
    while gen <= ngenerations
        % crossover
        pop = pop(randperm(popsize),:);  % randomize mating pool
        for i = 1 : popsize/2
           parent1index = i;
           parent2index = popsize - (i-1);
           [pop(parent1index,:),pop(parent2index,:)] = crossover(pop(parent1index,:),pop(parent2index,:),pc, xmin, xmax, stringlength, fun);
        end
        
        % mutation
        for j = 1 : popsize 
            pop(j,:) = mutation(pop(j,:),pm, xmin, xmax, stringlength, fun);
        end
        
        % selection
        pop = roulette(pop, stringlength, popsize);
        gen = gen + 1;
    end
    [optimum, index] = max(pop(:,stringlength+2)); 
    optimizer = pop(index,stringlength+1);

end


function [pop] = initialise(popsize, stringlength, a, b, fun)
    % randomly determine initial population
    pop = round(rand(popsize, stringlength+2));
    % calculate value of x and store in stringlength+1 th column
    bin2dvector = 2.^(size(pop(:,1:stringlength),2)-1:-1:0);
    pop(:,stringlength+1) = ((pop(:,1:stringlength)*(bin2dvector)')*(b-a))/(2.^stringlength-1)+a;
    % calculate function value and store in stringlenth+2 th column
    pop(:, stringlength+2) = arrayfun(fun, pop(:, stringlength+1));
end

function [child1, child2] = crossover(parent1, parent2, pc, a, b, stringlength, fun)
    if (rand<pc)
        % determine site for single point crossover
        cpoint = round(rand*(stringlength-1))+1;
        % perform crossover
        child1 = [parent1(:,1:cpoint) parent2(:,cpoint+1:stringlength)];
        child2 = [parent2(:,1:cpoint) parent1(:,cpoint+1:stringlength)];
        % update value of x after crossover in stringlength+1 th column
        bin2dvectorc1 = 2.^(size(child1(:,1:stringlength),2)-1:-1:0);
        bin2dvectorc2 = 2.^(size(child2(:,1:stringlength),2)-1:-1:0);
        child1(:,stringlength+1) = ((child1(:,1:stringlength)*(bin2dvectorc1)')*(b-a))/(2.^stringlength-1)+a;
        child2(:,stringlength+1) = ((child2(:,1:stringlength)*(bin2dvectorc2)')*(b-a))/(2.^stringlength-1)+a;
        % update function value and store in stringlength+2 th column
        child1(:, stringlength+2) = arrayfun(fun, child1(:, stringlength+1));
        child2(:, stringlength+2) = arrayfun(fun, child2(:, stringlength+1));
    else
        child1 = parent1;
        child2 = parent2;
    end
end

function [child] = mutation(parent, pm, a, b, stringlength, fun)
    if (rand < pm)
        % determine mutation site
        mpoint = round(rand*(stringlength-1))+1;
        child = parent;
        % flip the bit at mutation point
        child(mpoint) = abs(parent(mpoint)-1);
        % update value of x and store in stringlength+1 th column
        bin2dvectorc = 2.^(size(child(:,1:stringlength),2)-1:-1:0);
        child(:, stringlength+1) = ((child(:,1:stringlength)*(bin2dvectorc)')*(b-a))/(2.^stringlength-1)+a;
        % update function value and store in stringlength+2 th column
        child(:, stringlength+2) = arrayfun(fun, child(:, stringlength+1));
    else
        child = parent;
    end
end

function [newpop] = roulette(oldpop, stringlength, popsize)
    % calculate total fitness
    totalfit=sum(oldpop(:,stringlength+2));
    % calculate relative fitness
    prob=oldpop(:,stringlength+2) / totalfit;
    % sort by cumulative fitness
    prob = cumsum(prob);
    % generate probabilites for roulette wheel selection
    rns=sort(rand(popsize,1));
    fitin=1;    % index into old population
    newin=1;    % index into new population
    % select new population, keeping population size constant
    newpop = zeros(popsize, stringlength+2);
    while newin<=popsize 
        if (rns(newin)<prob(fitin))
            newpop(newin,:)=oldpop(fitin,:);
            newin = newin + 1;
        else
            fitin=fitin+1;
        end
    end
end


