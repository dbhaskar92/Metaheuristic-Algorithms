%
% Shifted Rotated Rastrigin's Function
% (xmin, xmax) = (-5, 5)
% error tolerance: 1
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = RastriginRotated(dim, seed)

	xmin = -5;
	xmax = 5;
	err = 1;

	% x*
	[fopt_loc, seed] = generate_random_vector(dim, seed, xmin, xmax);
	
    % f(x*) (interval randomly chosen)
	[value, seed] = r4_uni(seed);
	fopt = -211 + 422*value;
	
	% Random linear transformation matrix, condition number = 2
    cn = 2;
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

	function [ y ] = RastriginRotated_Def(x)

		z = (x - fopt_loc)*T;
		
		y = fopt;
		
		for i = 1 : dim
		    y = y + z(i)^2 - 10*cos(2*pi*z(i)) + 10;
		end
		
	end
	
	func = @RastriginRotated_Def;
	
end
