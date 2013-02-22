%% Steup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 100;
nNodes = 4;

totDiff = zeros(nTrials, 1);
trueLogZGaps = zeros(nTrials, 1);
lbpLogZGaps = zeros(nTrials, 1);
epsilon = 1e-6;

%% Loop it
for t = 1:nTrials
    %% Set up the problem    
    theta = -2*rand(nNodes, 1);    
    W = rand(nNodes, nNodes);
    W(1:nNodes+1:nNodes*nNodes) = 0;   
    W = .5 * (W + W');
    
    W = sparse(W);

    %% Run the algorithms
    [trueLogZ, ~, trueOneMarg, trueTwoMarg] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');   
    [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg]    = solveDAI(theta, W, 'BP', '[tol=1e-9,logdomain=0,updates=SEQRND]');    
    
    [CHECKlogZ, CHECKoneMarg, CHECKtwoMarg, CHECKmisc] = solveBethe(theta, W, epsilon);
    %[logZ, oneMarg, twoMarg, misc] = solveBethe(theta, W, epsilon);
    %intervalSz = getIntervalSz(misc.A, misc.B, W, epsilon);
    opts = struct();
    [logZ, oneMarg, twoMarg, misc] = BetheApprox_mex(theta, W, epsilon, opts);
    assertElementsAlmostEqual(logZ, CHECKlogZ);    
    assertElementsAlmostEqual(oneMarg, CHECKoneMarg);
    assertElementsAlmostEqual(twoMarg, CHECKtwoMarg);
    assertMiscEqual(misc, CHECKmisc);
    
    % The BBP bounds were only approximately calculated anyways
    %boundThresh = 0.002;
    %assert(all(betheMisc.A <= boundThresh + trueOneMarg & ...
    %           trueOneMarg - boundThresh <= (1 - betheMisc.B)), 'Bound violated!');
               
    %betheMisc.e    
    trueGap = trueLogZ - logZ;
    assert(trueGap > 0, 'Approximate Bethe logZ did not lower bound true logZ');
    lbpGap  = lbpLogZ - logZ;
    if lbpGap < 0
        warning('LBP logZ - Approximate Bethe logZ gap was negative');
    end
    
    if abs(lbpGap) > epsilon
        fn = ['questionable_case ' datestr(now, 0);];
        warning(['LBP logZ - Approximate Bethe logZ gap exceeded epsilon. ' ...
                 'Either LBP converged to a non-optimal stationary point ' ...
                 'or the algorithm has a bug. Example saved as in' fn]);                 
        save(fn);
    end          

    oneMargErr = mean(abs(trueOneMarg - oneMarg));
    disp(['Optimization finished with epsilon = ' num2str(epsilon)]);
    disp(['True logZ Gap = ' num2str(trueGap)]);
    disp(['LBP logZ Gap = ' num2str(lbpGap)]);
    disp(['Average error of one-marginals: ' num2str(oneMargErr) ]);    
    disp(['Average error of two-marginals: ' num2str(mean(abs(trueTwoMarg(:) - twoMarg(:)))) ]);    
        
    %% Record information
    trueLogZGaps(t)  = trueGap;
    lbpLogZGaps(t) = lbpGap;
    totDiff(t) = oneMargErr;
end

epsilonTxt = [' for \epsilon = ' num2str(epsilon)];

figure;
hist(totDiff);
title(['Average error of one-marginals' epsilonTxt]);

figure;
hist(trueLogZGaps);
title(['True logZ - Approximate Bethe logZ ' epsilonTxt]);

figure;
hist(lbpLogZGaps);
title(['LBP logZ - Approximate Bethe logZ ' epsilonTxt]);
