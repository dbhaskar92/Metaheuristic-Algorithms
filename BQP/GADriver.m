%
% Author: Dhananjay Bhaskar
% Script to solve UBQP using GA 
%

% params
problemfile = 'm7bqp500.txt';   % file containing matrix Q
timelimit = 120;                % time limit for each trial 
numtrials = 1;                  % number of trials
localsearch = 1;                % 0 : disable, 1 : enable

% initialize
tElapsed = zeros(1, numtrials);
optimum = zeros(1, numtrials);
gen = zeros(1, numtrials);

for trials = 1 : numtrials
    tStart = tic;
    [optimum(trials), gen(trials)] = GA(problemfile, timelimit, localsearch);
    tElapsed(trials) = toc(tStart);
end

% results
averageGen = mean(gen)
bestOpt = max(optimum)
averageOpt = mean(optimum)
averageTime = mean(tElapsed)