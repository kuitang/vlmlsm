%% Steup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 100;
nNodes = 6;

totDiff = zeros(nTrials, 1);
trueLogZGaps = zeros(nTrials, 1);
betheLogZGaps = zeros(nTrials, 1);
epsilon = 1e-4;

%% Loop it
for t = 1:nTrials
    %% Set up the problem    
    theta = -2*rand(nNodes, 1);    
    W = rand(nNodes, nNodes);
    W(1:nNodes+1:nNodes*nNodes) = 0;   
    W = .5 * (W + W');

    %% Run the algorithms
    [trueLogZ, joint, trueOneMarg, trueTwoMarg] = solveJTree(theta, W);   
    %betheLogZ = getBetheLogZ(theta, W, trueOneMarg, trueTwoMarg);
    [logZ, oneMarg, twoMarg, betheMisc] = solveBethe(theta, W, epsilon);
    
    % The BBP bounds were only approximately calculated anyways
    %boundThresh = 0.002;
    %assert(all(betheMisc.A <= boundThresh + trueOneMarg & ...
    %           trueOneMarg - boundThresh <= (1 - betheMisc.B)), 'Bound violated!');
               
    %betheMisc.e
    %betheGap = betheLogZ - logZ;
    trueGap  = trueLogZ  - logZ;
    assert(trueGap > 0, 'Approximate Bethe logZ did not lower bound true logZ');
    %assert(betheGap > epsilon, 'Approximate Bethe logZ not within epsilon of true Bethe logZ');    

    oneMargErr = mean(abs(trueOneMarg - oneMarg));
    disp(['Optimization finished with epsilon = ' num2str(epsilon)]);
    disp(['True logZ Gap = ' num2str(trueGap)]);
    disp(['Average error of one-marginals: ' num2str(oneMargErr) ]);    
    disp(['Average error of two-marginals: ' num2str(mean(abs(trueTwoMarg(:) - twoMarg(:)))) ]);    
        
    %% Record information
    trueLogZGaps(t)  = trueGap;
    betheLogZGaps(t) = betheGap;
    totDiffs(t) = oneMargErr;
end

epsilonTxt = [' for \epsilon = ' num2str(epsilon)];

figure;
hist(totDiffs);
title(['Average error of one-marginals' epsilonTxt]);

figure;
hist(trueLogZGaps);
title(['True logZ - Approximate Bethe logZ ' epsilonTxt]);

% figure;
% hist(betheLogZGaps);
% title(['Bethe logZ - Approximate Bethe logZ ' epsilonTxt]);
