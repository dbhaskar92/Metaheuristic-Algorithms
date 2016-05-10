%
% Shifted Rosenbrocks's Function
% (xmin, xmax) = (-30, 30)
% error tolerance: 1
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = Rosenbrock(dim, seed)

	xmin = -30;
	xmax = 30;
	err = 1;

	% x*
	[fopt_loc, seed] = generate_random_vector(dim, seed, xmin, xmax);
	
    % f(x*) (interval randomly chosen)
	[value, seed] = r4_uni(seed);
	fopt = -390 + 780*value;

	function [ y ] = Rosenbrock_Def(x)
		
		z = x - fopt_loc + 1;
		
		y = fopt;
		
		for i = 1 : dim - 1
		    y = y + 100*((z(i)^2 - z(i+1))^2) + (z(i) - 1)^2;
		end
		
	end

	func = @Rosenbrock_Def;

end
