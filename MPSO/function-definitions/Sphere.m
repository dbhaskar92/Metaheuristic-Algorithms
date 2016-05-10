%
% Shifted Sphere Function
% (xmin, xmax) = (-100, 100)
% error tolerance: 1e-3
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = Sphere(dim, seed)

	xmin = -100;
	xmax = 100;
	err = 0.001;

	% x*
	[fopt_loc, seed] = generate_random_vector(dim, seed, xmin, xmax);
	
	% f(x*) (interval randomly chosen)
	[value, seed] = r4_uni(seed);
	fopt = -450 + 900*value;
	
	function [ y ] = Sphere_Def(x)
		
		z = x - fopt_loc;
		
		y = fopt;
		
		for i = 1 : dim
		    y = y + z(i)^2;
		end
		
	end
	
	func = @Sphere_Def;
	
end
