%
% Shifted Rotated Griewank's Function
% (xmin, xmax) = (-600, 600)
% error tolerance: 1
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = Griewank(dim, seed)

	xmin = -600;
	xmax = 600;
	err = 1;

	% x*
	[fopt_loc, seed] = generate_random_vector(dim, seed, xmin, xmax);
	
	% f(x*) (interval randomly chosen)
	[value, seed] = r4_uni(seed);
	fopt = -180 + 360*value;
	
	% Random linear transformation matrix, condition number 3
	cn = 3;
    T_cn = 0;
    T = zeros(dim, dim);
    while T_cn ~= cn
        [R, seed] = generate_random_matrix(dim, dim, seed);
        [U, S, V] = svd(R);
        S(S ~= 0) = linspace(cn, 1, min(dim, dim));
        T = U*S*V';
        T_cn = cond(T);
    end
    assert(cond(T)==cn, 'Failed to generate matrix with required condition number');
    
	function [ y ] = Griewank_Def(x)
	
		z = (x - fopt_loc)*T;

		y = fopt;
		
		for i = 1 : dim
		    y = y + (z(i)^2)/4000;
		end
		prod = 1;
		for i = 1 : dim
		    prod = prod*cos(z(i)/sqrt(i));
		end
		y = y - prod + 1;
		
	end
	
	func = @Griewank_Def;
	
end
