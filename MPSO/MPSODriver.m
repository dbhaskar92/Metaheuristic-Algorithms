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
    
    [f_handle, fopt, fopt_loc, xmin, xmax, err, seed] = feval(fun, 2, seed);
    
    Plot3DObjective(f_handle, char(func_names(i)), fopt, fopt_loc, xmin, xmax, func_plot_stepsize(i));
    
end