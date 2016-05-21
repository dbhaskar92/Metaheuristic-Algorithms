%
% Description: Driver script to plot benchmark functions and perform MPSO runs 
% in high dimensional search space
%
% Author: Dhananjay Bhaskar <dbhaskar92@gmail>
%
% Comments: We use ZIGGURAT PRNG (original C code written by
% Marsaglia and Tsang, translated to MATLAB by John Burkardt) to generate
% random numbers from uniform and normal distributions
%

addpath('./random/');
addpath('./function-definitions/');
seed = 123;

func_names = {'Ackley' 'Elliptic' 'Griewank' 'Rastrigin' 'RastriginRotated' 'Rosenbrock' 'Schwefel' 'Sphere' 'Weierstrass'};
func_plot_stepsize = [0.5 1 10 0.05 0.05 0.1 2 2 0.01];

% Plot 3-D maps of all benchmark functions
for i = 1 : numel(func_names)
    
    fun = str2func(char(func_names(i)));
    
    [f_handle, fopt, fopt_loc, xmin, xmax, ~, seed] = feval(fun, 2, seed);
    
    Plot3DObjective(f_handle, char(func_names(i)), fopt, fopt_loc, xmin, xmax, func_plot_stepsize(i));
    
end

diary('output.txt')

% Trial Runs

N = 15;
dim = 5;
[f_handle, fopt, fopt_loc, xmin, xmax, err, seed] = feval('Ackley', dim, seed);

% RWDE params
% Scheme:
%   0: no local search
%   1: local search on global best (gbest) candidate solution
%   2: local search on gbest and pbest of randomly selected particles
scheme = 2;
freq = 10;
eps = 0.1;
tmax = 5;
lbd = 1.0;

% RWMPSOg
[opt, res, fe, seed] = MPSOn(N, lbd, eps, tmax, freq, xmin, xmax, scheme, f_handle, dim, err, fopt, seed);
fprintf('\nTrial Run of RWMPSOg for Ackley Function (10000 iterations): \n');
fprintf('EXPECTED Optimum: %f Location: %s \n', fopt, mat2str(fopt_loc));
fprintf('RESULTS Optimum: %f Location: %s FEV: %d \n', opt, mat2str(res), fe);

% RWMPSOl
[opt, res, fe, seed] = MPSOln(N, lbd, eps, tmax, freq, xmin, xmax, scheme, f_handle, dim, err, fopt, seed);
fprintf('\nTrial Run of RWMPSOl for Ackley Function (10000 iterations): \n');
fprintf('EXPECTED Optimum: %f Location: %s \n', fopt, mat2str(fopt_loc));
fprintf('RESULTS Optimum: %f Location: %s FEV: %d \n', opt, mat2str(res), fe);

% Empirical Experiments

dim = 10;
pop_size = [15 30 60];
num_trials = 10;

for i = 1 : numel(func_names)
    
    fun = str2func(char(func_names(i)));
    
    for j = 1 : numel(pop_size)
        
        N = pop_size(j);
        
        % RWMPSOg
        optimum = zeros(num_trials, 1);
        fev = zeros(num_trials, 1);
        fail = 0;

        for k = 1 : num_trials
            
            [fh, fopt, floc, xmin, xmax, err, seed] = feval(fun, dim, seed);
            [optimum(k,1), ~, fev(k,1), seed] = MPSOn(N, lbd, eps, tmax, freq, xmin, xmax, scheme, fh, dim, err, fopt, seed);
            
            if (abs(optimum(k,1) - fopt) > err)
                fev(k,1) = 0;
                fail = fail + 1;
            end
            
        end

        fev(fev==0) = [];
        stats = [min(fev(:,1)) mean(fev(:,1)) max(fev(:,1)) std(fev(:,1))];
        fprintf('\nRWMPSOg Empirical results for %s \n', char(func_names(i)));
        fprintf('SS: %d Min/Mean/Max/Std: %s Suc: %d \n', N, mat2str(stats), num_trials-fail);
        
        % RWMPSOl
        optimum = zeros(num_trials, 1);
        fev = zeros(num_trials, 1);
        fail = 0;

        for k = 1 : num_trials
            
            [fh, fopt, floc, xmin, xmax, err, seed] = feval(fun, dim, seed);
            [optimum(k,1), ~, fev(k,1), seed] = MPSOln(N, lbd, eps, tmax, freq, xmin, xmax, scheme, fh, dim, err, fopt, seed);
            
            if (abs(optimum(k,1) - fopt) > err)
                fev(k,1) = 0;
                fail = fail + 1;
            end
            
        end

        fev(fev==0) = [];
        stats = [min(fev(:,1)) mean(fev(:,1)) max(fev(:,1)) std(fev(:,1))];
        fprintf('\nRWMPSOl Empirical results for %s \n', char(func_names(i)));
        fprintf('SS: %d Min/Mean/Max/Std: %s Suc: %d \n', N, mat2str(stats), num_trials-fail);
        
    end
    
end
