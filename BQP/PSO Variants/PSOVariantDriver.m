%
% Driver to test various PSO variants on problems taken from ORLIB
%
% Author: Dhananjay Bhaskar
%
% Last modified: Saturday, Oct 26, 2013
%

% ORLIB matrix
problemfile = 'm6bqp50.txt';

% number of generations
ngen = [500, 1000, 1500, 2000];

% population size
popsize = [30, 60, 120];

for p = 1 : numel(ngen)
    for q = 1 : numel(popsize)

        [~, optimum1, timetot1] = MBPSO(problemfile, ngen(p), popsize(q));
        
        fprintf('\n MBPSO - Ngen : %i Popsize : %i Optimum : %i Time : %d', ngen(p), popsize(q), optimum1, timetot1);

        [~, optimum2, timetot2] = KBPSO(problemfile, ngen(p), popsize(q));
        
        fprintf('\n KBPSO - Ngen : %i Popsize : %i Optimum : %i Time : %d \n', ngen(p), popsize(q), optimum2, timetot2);
        
    end
end