%
% Ackley's Function
% (xmin, xmax) = (-32, 32)
% error tolerance: 1e-3
%

function [func, fopt, fopt_loc, xmin, xmax, err, seed] = Ackley(dim, seed)

	xmin = -32;
	xmax = 32;
	err = 0.001;
	
	fopt = 0;
	fopt_loc = zeros(1,dim);

	function [ y ] = Ackley_Def(x)
		
		arg1 = 0;
		for i = 1 : dim
		    arg1 = arg1 + x(i)^2;
		end
		arg1 = arg1/dim;
		arg1 = -0.2*sqrt(arg1);
		
		arg2 = 0;
		for i = 1 : dim
		    arg2 = arg2 + cos(2*pi*x(i));
		end
		arg2 = arg2/dim;
		
		y = 20 - 20*exp(arg1) - exp(arg2) + exp(1);
	end
		
	func = @Ackley_Def;
	
end
