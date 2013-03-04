%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 10;

nNodes = 1000;
runJT = true;
testTrees = true;
checkMatlab = true;

epsilon = 1e-2;

%% Zero-out some vectors... does it even matter?
% totDiff(nTrials) = 0;
% lbpTotDiff(nTrials) = 0;
% 
% trueLogZGaps(nTrials) = 0;
% lbpLogZGaps(nTrials) = 0;
% 
% betheTimes(nTrials) = 0;
% bkEdges(nTrials) = 0;
% maxFlowTimes(nTrials) = 0;
% lbpTimes(nTrials) = 0;
% JTTimes(nTrials) = 0;

%% Loop it
for t = 1:nTrials
    %% Set up the problem
    if testTrees
        [theta, W] = makeProblem(nNodes, -2, 1);
        T = randTree(nNodes);
        W = W .* T;
    else
        [theta, W] = makeProblem(nNodes, -0.5, 0);
    end
        
    %% Run the algorithms
    opts = struct();
    
    tic;
    [logZ, oneMarg, twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, opts);        
    betheTimes(t) = toc;
    
    bkEdges(t) = misc.nBKEdges;
    maxFlowTimes(t)  = misc.BKMaxFlowTime;
    
    times = [misc.makeMinSumTime misc.BKConstructTime misc.BKMaxFlowTime] ./ misc.mexTotTime;
    fprintf(1, 'Make minsum = %g, BK construction = %g, Max flow = %g, Overhead = %g\n', ...
            times(1), times(2), times(3), 1 - sum(times));
    fprintf(1, 'MEX-call overhead: %g\n', 1 - misc.mexTotTime / betheTimes(t));    
        
    misc
    
    if runJT        
        [trueLogZ, ~, trueOneMarg, trueTwoMarg, JTTimes(t)] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');           
    end
    
    tic;
    [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, lbpTimes(t)] = solveDAI(theta, W, 'BP', '[tol=1e-9,logdomain=0,updates=SEQRND]');        
    
    if checkMatlab
        %% Check against MATLAB (cell execution)
        [CHECKlogZ, CHECKoneMarg, CHECKtwoMarg, CHECKmisc] = solveBethe(theta, W, epsilon);
        [logZ, oneMarg, twoMarg, misc] = BetheApprox_debug_mex(theta, W, epsilon, opts);
        intervalSz = getIntervalSz(misc.A, misc.B, W, epsilon);

        assertElementsAlmostEqual(logZ, CHECKlogZ);    
        assertElementsAlmostEqual(oneMarg, CHECKoneMarg);
        assertElementsAlmostEqual(twoMarg, CHECKtwoMarg);
        assertMiscEqual(misc, CHECKmisc);
    end
    
    lbpGap  = lbpLogZ - logZ;
    if lbpGap < 0
        warning('LBP logZ - Approximate Bethe logZ gap was negative');
    end

    lbpOneMargErr = mean(abs(lbpOneMarg - oneMarg));
    disp(['LBP logZ Gap = ' num2str(lbpGap)]);               
    disp(['Average lbp-Bethe deviation of one-marginals: ' num2str(lbpOneMargErr) ]);
    disp(['Average lbp-Bethe deviation of two-marginals: ' num2str(mean(abs(lbpTwoMarg(:) - twoMarg(:)))) ]);    
    lbpLogZGaps(t) = lbpGap;
    lbpTotDiff(t) = lbpOneMargErr;
    
    if abs(lbpGap) > epsilon
        fn = ['questionable_case ' datestr(now, 0);];
        warning(['LBP logZ - Approximate Bethe logZ gap exceeded epsilon. ' ...
                 'Either LBP converged to a non-optimal stationary point ' ...
                 'or the algorithm has a bug. Example saved as in' fn]);                 
        save(fn);
    end          

    disp(['Optimization finished with epsilon = ' num2str(epsilon)]);
    
    if runJT
        trueGap = trueLogZ - logZ;
        assert(trueGap > 0, 'Approximate Bethe logZ did not lower bound true logZ');
        
        oneMargErr = mean(abs(trueOneMarg - oneMarg));
        
        disp(['True logZ Gap = ' num2str(trueGap)]);
        disp(['Average error of one-marginals: ' num2str(oneMargErr) ]);
        disp(['Average error of two-marginals: ' num2str(mean(abs(trueTwoMarg(:) - twoMarg(:)))) ]);
        trueLogZGaps(t)  = trueGap;
        totDiff(t) = oneMargErr;
    end

end

betheTimes
lbpTimes
if runJT
    JTTimes
end

figure;
[~, idxs] = sort(bkEdges);
plot(bkEdges(idxs), maxFlowTimes(idxs), bkEdges(idxs), betheTimes(idxs));
legend('Max Flow Time', 'Total Bethe Time');
title('Runtime vs # of BK Edges');

epsilonTxt = [' for \epsilon = ' num2str(epsilon)];

figure;
hist(double(bkEdges));
title(['Number of BK Edges' epsilonTxt]);

if runJT
    figure;
    hist(totDiff);
    title(['Average error of one-marginals' epsilonTxt]);
    
    figure;
    hist(trueLogZGaps);
    title(['True logZ - Approximate Bethe logZ ' epsilonTxt]);

    figure;
    hist(lbpLogZGaps);
    title(['LBP logZ - Approximate Bethe logZ ' epsilonTxt]);
end

