%
% Generate a random matrix of given dimensions
% Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
%

function [res, seed] = generate_random_matrix(nrows, ncols, seed)

	res = zeros(nrows, ncols);
	
	[kn, fn, wn] = r4_nor_setup();

	for i = 1 : nrows
	
		for j = 1 : ncols
		
			[value, seed] = r4_nor(seed, kn, fn, wn);
			
			res(i, j) = value;
		
		end
		
	end	

end
