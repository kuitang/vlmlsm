%% Steup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 100;
nNodes = 5;

totDiffs = zeros(nTrials, 1);
logZGaps = zeros(nTrials, 1);
epsilon = 1e-4;

%% Loop it
for t = 1:nTrials
    %% Set up the problem    
    theta = -2*rand(nNodes, 1);    
    W = rand(nNodes, nNodes);
    W(1:nNodes+1:nNodes*nNodes) = 0;   
    W = .5 * (W + W');

    %% Run the algorithms
    [trueLogZ, joint, trueOneMarginals, trueTwoMarginals] = solveJTree(theta, W);   
    [logZ, oneMarginals, twoMarginals] = solveBethe(theta, W, epsilon);
    trueOneMarginals
    
    totDiffs(t) = sum(abs(trueOneMarginals - oneMarginals));
    gap = trueLogZ - logZ;
    assert(gap > 0,       'Uh oh, Bethe logZ did not lower bound true logZ');
    %assert(gap < epsilon, 'Uh oh, Bethe logZ is not within epsislon of true logZ');
    logZGaps(t) = gap;
    
    disp(['Optimization finished with epsilon = ' num2str(epsilon)]);
    disp(['One-norm of one-marginals: ' num2str(sum(abs(trueOneMarginals - oneMarginals))) ]);
    disp(['Two-norm of one-marginals: ' num2str(norm(trueOneMarginals - oneMarginals)) ]);
    disp(['One-norm of two-marginals: ' num2str(sum(abs(trueTwoMarginals(:) - twoMarginals(:)))) ]);
    disp(['Two-norm of two-marginals (lineralized): ' num2str(norm(trueTwoMarginals(:) - twoMarginals(:))) ]);
end

figure;
hist(totDiffs);
title(['1-norm differences of brute force and computed marginals for \epsilon = ' num2str(epsilon)]);

figure;
hist(logZGaps);
title(['True logZ - Bethe logZ for \epsilon = ' num2str(epsilon)]);
