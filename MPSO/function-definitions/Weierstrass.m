%
% Shifted Rotated Weierstrass Function
% (xmin, xmax) = (-0.5, 0.5)
% error tolerance: 1
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = Weierstrass(dim, seed)

	xmin = -0.5;
	xmax = 0.5;
	err = 1;

	% x*
	[fopt_loc, seed] = generate_random_vector(dim, seed, xmin, xmax);
	
    % f(x*) (interval randomly chosen)
	[value, seed] = r4_uni(seed);
	fopt = -90 + 180*value;
	
	% Random linear transformation matrix, condition number = 5
    cn = 5;
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
	
	function [ y ] = Weierstrass_Def(x)
	
		a = 0.5;
		b = 3;
		kmax = 20;
		
		z = (x - fopt_loc)*T;
		
		y = fopt;
		
		for i = 1 : dim
		    sum1 = 0;
		    for k = 0 : kmax
		        sum1 = sum1 + (a^k)*cos(2*pi*(b^k)*(z(i)+0.5));
		    end
		    y = y + sum1;
		end
		sum2 = 0;
		for k = 0 : kmax
		    sum2 = sum2 + (a^k)*cos(2*pi*(b^k)*0.5);
		end
		y = y - dim*sum2;
		
	end
	
	func = @Weierstrass_Def;
	
end
