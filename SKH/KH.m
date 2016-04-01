% 
% Implementation of KH (Krill herd algorithm)
%
% Author:
%   Dhananjay Bhaskar
%
% References:
%   Krill herd: A new bio-inspired optimization algorithm
%   Gandomi, Alavi
%
% Last modified:
%   Sunday, Aug 18, 2013
%

function [minima, minimizer] = KH(xmin, xmax, fun, dim)

    % params
    N = 50;         % population size
    MaxGen = 50;    % maximum number of generations
    Vf = 0.02;      % foraging speed
    Dmax = 0.005;   % maximum diffusion speed
    Nmax = 0.01;    % maximum induced speed
    Keep = 15;      % maximum number of acceptances
    elitism = 0;    % 0: disable 1: enable
    
    % initialization
    t = 0;                          % generation counter
    Ct = 0.5;                       % constant lies in [0,2]
    deltaT = Ct*(xmax-xmin)*dim;    % step size 
    Nim = zeros(N, dim);            % induced motion
    F = zeros(N, dim);              % foraging motion
   
    x = xmin + (xmax - xmin).*rand(N,dim+1);    % initialize position
    
    % fitness calculation
    for i = 1 : N
        x(i,dim+1) = feval(fun, x(i,1:dim));
    end
    
    % keep track of best positions visited by each krill
    indbest = zeros(N,dim+2);
    indbest(:,1:dim) = x(:,1:dim);
    for i = 1 : N
        indbest(i,dim+1) = x(i, dim+1);  % store curr fitness in col 1
        indbest(i,dim+2) = x(i, dim+1);  % store best fitness in col 2
    end
    
    while (t < MaxGen)
        
            % sort by fitness
            x = sortrows(x,dim+1);
            indbest = sortrows(indbest, dim+1);
            
            if (elitism == 1)
                % store the KEEP best krill
                bestKrills = zeros(Keep, dim+1);
                for i = 1 : Keep
                    bestKrills(i,:) = x(i,:);
                end
            end
            
            % krill herd update
            xnew = zeros(N,dim);
            
            for i = 1 : N
                
                % calculate Xfood (mean)
                Xfood = zeros(1, dim);
                for j = 1 : N
                    Xfood = Xfood + x(j,1:dim).*(1/(x(j,dim+1)));
                end
                denominator = 0;
                for j = 1 : N
                    denominator = denominator + (1/(x(j,dim+1)));
                end
                Xfood = Xfood ./ denominator;
                XfoodFit = feval(fun, Xfood);
                
                beta = calculateBetaFood(x, dim, Xfood, XfoodFit, i, t, MaxGen) + calculateBetaBest(x, i, dim, indbest); 
                omegaForaging = (45 - 0.8*t) / 50;
                F(i,:) = Vf.*beta + omegaForaging.*F(i,:);          % foraging motion
                
                omegaInduced = (45 - 0.8*t) / 50;
                alpha = calculateAlphaLoc(x, dim, i, N) + calculateAlphaTarget(x, indbest, dim, i, t, MaxGen);
                Nim(i,:) = Nmax.*alpha + omegaInduced.*Nim(i,:);    % induced motion
                
                delta = -1 + (1+1).*rand(1,dim);                    % random directional vector
                D = Dmax.*delta;                                    % physical diffusion
                
                xnew(i,:) = x(i,1:dim) + deltaT.*(F(i,:) + Nim(i,:) + D);
                
                x(i, 1:dim) = xnew(i,:);
             
                x(i,dim+1) = feval(fun, x(i,1:dim));                % update fitness
               
            end
            
            % update individual best positions
            for i = 1 : N
                indbest(i,dim+1) = x(i, dim+1);
                if (x(i, dim+1) < indbest(i, dim+2))
                    indbest(i,1:dim) = x(i,1:dim);
                    indbest(i,dim+2) = x(i,dim+1);
                end
            end
            
            % elitest strategy
            if (elitism == 1)
                x = sortrows(x, -1*(dim+1));                            % sort in descensding order
                indbest = sortrows(indbest, -(dim+1));
                for i = 1 : Keep
                    x(i,:) = bestKrills(i,:);                           % replace Keep worst krills
                    indbest(i,1:dim) = x(i,1:dim);
                    indbest(i,dim+1) = x(i,dim+1);
                    indbest(i,dim+2) = x(i,dim+1);
                end
            end
            
            % update
            t = t + 1;
    end
    
    [minima, index] = min(x(:,dim+1));
    minimizer = x(index, 1:dim);
    
end

function [unitX] = Xhat(x, dim, p, q)
        numerator = x(q,1:dim) - x(p, 1:dim);
        denominator = norm(numerator,2) + 0.0001;
        unitX = numerator ./ denominator;
end

function [unitK] = Khat(x, dim, p, q)
        numerator = x(p,dim+1) - x(q,dim+1);
        [Kbest, ~] = min(x(:,dim+1));
        [Kworst, ~] = max(x(:,dim+1));
        denominator = Kworst - Kbest;
        unitK = numerator / denominator;
end

function [alphaloc] = calculateAlphaLoc(x, dim, i, N)
        alphaloc = zeros(1, dim);
        
        % calculate sensing distance
        dsensing = 0;
        for j = 1 : N
            dsensing = dsensing + norm(x(i,1:dim)-x(j,1:dim),2);
        end
        dsensing = dsensing / (5*N);
        
        for j = 1 : N
            neighbourdist = norm(x(i,1:dim)-x(j,1:dim),2);
            if (neighbourdist < dsensing)
                alphaloc = alphaloc + Khat(x, dim, i, j).*Xhat(x, dim, i, j);
            end
        end
end

function [alphatarget] = calculateAlphaTarget(x, indbest, dim, i, t, MaxGen)
        
        Cbest = rand(1) + ((t+1)/MaxGen);
        Cbest = Cbest * 2;
        
        [~, best] = min(indbest(:,dim+2));
        alphatarget = Khat(x, dim, i, best).*Xhat(x, dim, i, best);
        
        alphatarget = Cbest.*alphatarget;
end

function [betafood] = calculateBetaFood(x, dim, Xfood, XfoodFit, i, t, MaxGen)
        
        Cfood = 1 - ((t+1)/MaxGen);
        Cfood = Cfood * 2;
        
        % calculate X i,food
        numer = Xfood - x(i, 1:dim);
        denom = norm(numer,2) + 0.0001;
        Xunit = numer ./ denom;
        
        % calculate K i,food
        numer1 = x(i,dim+1) - XfoodFit;
        [Kbest, ~] = min(x(:,dim+1));
        [Kworst, ~] = max(x(:,dim+1));
        denom1 = Kworst - Kbest;
        Kunit = numer1 / denom1;
        
        betafood = Kunit.*Xunit;
        betafood = Cfood.*betafood;
end

function [betabest] = calculateBetaBest(x, i, dim, indbest)
        
        % calculate X i,ibest
        numer = indbest(i,1:dim) - x(i, 1:dim);
        denom = norm(numer,2) + 0.0001;
        Xunit = numer ./ denom;
        
        % calculate K i,ibest
        numer1 = x(i,dim+1) - indbest(i,dim+2);
        [Kbest, ~] = min(x(:,dim+1));
        [Kworst, ~] = max(x(:,dim+1));
        denom1 = Kworst - Kbest;
        Kunit = numer1 / denom1;
        
        betabest = Kunit.*Xunit;
end

