%
% Generate a random vector of given dimension with each component in (cmin, cmax) 
% Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
%

function [res, seed] = generate_random_vector(dim, seed, cmin, cmax)
	
	res = zeros(1, dim);

	for i = 1 : dim
	
		[value, seed] = r4_uni(seed);
		
		res(i) = cmin + (cmax - cmin)*value;
	
	end

end
