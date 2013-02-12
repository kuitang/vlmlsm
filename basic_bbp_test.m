%% Steup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 100;
nNodes = 3;

totDiffs = zeros(nTrials, 1);
logZGaps = zeros(nTrials, 1);
epsilon = 1e-6;

%% Loop it
for t = 1:nTrials
    %% Set up the problem    
    theta = -2*rand(nNodes, 1);    
    W = rand(nNodes, nNodes);
    W(1:nNodes+1:nNodes*nNodes) = 0;   
    W = .5 * (W + W');

    %% Run the algorithms
    [trueLogZ, joint, trueOneMarginals, trueTwoMarginals] = solveJTree(theta, W);   
    [logZ, oneMarginals, twoMarginals, betheMisc] = solveBethe(theta, W, epsilon);
    
    % The BBP bounds were only approximately calculated anyways
    boundThresh = 0.002;
    assert(all(betheMisc.A <= boundThresh + trueOneMarginals & ...
               trueOneMarginals - boundThresh <= (1 - betheMisc.B)), 'Bound violated!');
    
    gap = trueLogZ - logZ;
    assert(gap > 0,       'Uh oh, Bethe logZ did not lower bound true logZ');
    if gap > epsilon
        warning('Uh oh, Bethe logZ is not within epsislon of true logZ');
    end


    disp(['Optimization finished with epsilon = ' num2str(epsilon)]);
    disp(['logZ Gap = ' num2str(gap)]);
    disp(['Average error of one-marginals: ' num2str(mean(abs(trueOneMarginals - oneMarginals))) ]);
    %disp(['Two-norm of one-marginals: ' num2str(norm(trueOneMarginals - oneMarginals)) ]);
    disp(['Average error of two-marginals: ' num2str(mean(abs(trueTwoMarginals(:) - twoMarginals(:)))) ]);
    %disp(['Two-norm of two-marginals (lineralized): ' num2str(norm(trueTwoMarginals(:) - twoMarginals(:))) ]);
        
    %% Record information
    logZGaps(t) = gap;
    totDiffs(t) = sum(abs(trueOneMarginals - oneMarginals));
end

figure;
hist(totDiffs);
title(['Average error of one-marginals for \epsilon = ' num2str(epsilon)]);

figure;
hist(logZGaps);
title(['True logZ - Bethe logZ for \epsilon = ' num2str(epsilon)]);
