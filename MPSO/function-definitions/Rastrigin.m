%
% Shifted Rastrigin's Function
% (xmin, xmax) = (-5.12, 5.12)
% error tolerance: 1
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = Rastrigin(dim, seed)

	xmin = -5.12;
	xmax = 5.12;
	err = 1;

	% x*
	[fopt_loc, seed] = generate_random_vector(dim, seed, xmin, xmax);
	
	% f(x*) (interval randomly chosen)
    [value, seed] = r4_uni(seed);
	fopt = -330 + 660*value;

	function [ y ] = Rastrigin_Def(x)
		
		z = x - fopt_loc;
		
		y = fopt;
		
		for i = 1 : dim
		    y = y + z(i)^2 - 10*cos(2*pi*z(i)) + 10;
		end
		
	end
	
	func = @Rastrigin_Def;
	
end
