%
% Shifted Rotated High Conditioned Elliptic Function
% (xmin, xmax) = (-100, 100)
% error tolerance: 500
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = Elliptic(dim, seed)

	xmin = -100;
	xmax = 100;
	err = 500;

	% x*
	[fopt_loc, seed] = generate_random_vector(dim, seed, xmin, xmax);
	
	% f(x*) (interval randomly chosen)
	[value, seed] = r4_uni(seed);
	fopt = -145 + 290*value;
	
	% random orthogonal matrix
	[R, seed] = generate_random_matrix(dim, dim, seed);
	R = orth(R);

	function [ y ] = Elliptic_Def(x)
		
		z = (x - fopt_loc)*R;
		
		y = fopt;
		
		for i = 1 : dim
		    y = y + (z(i)^2)*(10^6)^((i-1)/(dim-1));
		end
		
	end
	
	func = @Elliptic_Def;
	
end
