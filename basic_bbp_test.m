%% Steup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 100;
nNodes = 12;

runJT = false;
checkMatlab = false;

totDiff = zeros(nTrials, 1);
trueLogZGaps = zeros(nTrials, 1);
lbpLogZGaps = zeros(nTrials, 1);
epsilon = 1e-3;

%% Loop it
for t = 1:nTrials
    %% Set up the problem    
    theta = -nNodes/2*rand(nNodes, 1);
    W = makeMonge(nNodes, nNodes);    
    W(1:nNodes+1:nNodes*nNodes) = 0; 
    W = .5 * (W + W');
    W = W ./ max(W(:));    
    
    % Scale    
    W = sparse(W);    

    %% Run the algorithms
    opts = struct();

    [logZ, oneMarg, twoMarg, misc] = BetheApprox_mex(theta, W, epsilon, opts);    
    
    if runJT
        [trueLogZ, ~, trueOneMarg, trueTwoMarg] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');   
    end
    
    [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg]    = solveDAI(theta, W, 'BP', '[tol=1e-6,logdomain=0,updates=SEQRND]');    
    
    if checkMatlab
        [CHECKlogZ, CHECKoneMarg, CHECKtwoMarg, CHECKmisc] = solveBethe(theta, W, epsilon);
        [logZ, oneMarg, twoMarg, misc] = solveBethe(theta, W, epsilon);
        intervalSz = getIntervalSz(misc.A, misc.B, W, epsilon);

        assertElementsAlmostEqual(logZ, CHECKlogZ);    
        assertElementsAlmostEqual(oneMarg, CHECKoneMarg);
        assertElementsAlmostEqual(twoMarg, CHECKtwoMarg);
        assertMiscEqual(misc, CHECKmisc);
    end
    
    % The BBP bounds were only approximately calculated anyways
    %boundThresh = 0.002;
    %assert(all(betheMisc.A <= boundThresh + trueOneMarg & ...
    %           trueOneMarg - boundThresh <= (1 - betheMisc.B)), 'Bound violated!');
                   
    
    if runJT
        trueGap = trueLogZ - logZ;
        assert(trueGap > 0, 'Approximate Bethe logZ did not lower bound true logZ');
    end
    
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

    disp(['Optimization finished with epsilon = ' num2str(epsilon)]);
    
    if runJT
        disp(['True logZ Gap = ' num2str(trueGap)]);
        disp(['Average error of one-marginals: ' num2str(oneMargErr) ]);
        disp(['Average error of two-marginals: ' num2str(mean(abs(trueTwoMarg(:) - twoMarg(:)))) ]);
        oneMargErr = mean(abs(trueOneMarg - oneMarg));
        trueLogZGaps(t)  = trueGap;
        totDiff(t) = oneMargErr;        
    end    
    
    disp(['LBP logZ Gap = ' num2str(lbpGap)]);               
    lbpLogZGaps(t) = lbpGap;
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
