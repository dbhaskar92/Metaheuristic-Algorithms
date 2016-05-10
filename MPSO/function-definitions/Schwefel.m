%
% Shifted Schwefel's Function
% (xmin, xmax) = (-100, 100)
% error tolerance: 1e-2
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = Schwefel(dim, seed)

	xmin = -100;
	xmax = 100;
	err = 0.01;

	% x*
	[fopt_loc, seed] = generate_random_vector(dim, seed, xmin, xmax);
	
    % f(x*) (interval randomly chosen)
	[value, seed] = r4_uni(seed);
	fopt = -310 + 620*value;

	function [ y ] = Schwefel_Def(x)
		
		z = x - fopt_loc;
		
		y = fopt;
		
		for i = 1 : dim
		    tmp = 0;
		    for j = 1 : i
		        tmp = tmp + z(j);
		    end
		    y = y + tmp^2;
		end
		
	end
	
	func = @Schwefel_Def;
	
end
