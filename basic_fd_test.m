% Basic first derivative bound test.

%% Setup
nTrials = 10;
nNodes  = 16;

%nNodes = 10000;
runJT = true;
%testTrees = true;
testTrees = false;
treeLoops = 200;
checkMatlab = false;

% Uniform bounds for weights
ta = -2;
tb = 2;
wa = -1;
wb = 1;

density = 1;

epsilon = 1e-2;

oneMarg(nNodes,nTrials,2) = 0;
trueOneMarg(nNodes,nTrials,2) = 0;
A(nNodes,nTrials,2) = 0;
B(nNodes,nTrials,2) = 0;
ABgap(nNodes,nTrials,2) = 0;
intervalSzs(nTrials,2) = 0;

% TODO: Test various sparsity patterns.

%% Zero-out some vectors... does it even matter?
totDiff(nTrials,2) = 0;
lbpTotDiff(nTrials,2) = 0;

trueLogZGaps(nTrials,2) = 0;
lbpLogZGaps(nTrials,2) = 0;

betheTimes(nTrials,2) = 0;
bkEdges(nTrials,2) = 0;
maxFlowTimes(nTrials,2) = 0;

lbpTimes(nTrials) = 0;
JTTimes(nTrials) = 0;

opts = struct('useMooij', false);
mooijOpts = struct('useMooij', true);

problems = cell(1, nTrials);

%% Loop it
for t = 1:nTrials
    %% Set up the problem
    if testTrees
        T = randTree(nNodes, nLoops);
        [theta, W] = makeUnifProblemNew(nNodes, T, ta, tb, wa, wb)
    else
        [theta, W] = makeUnifProblemNew(nNodes, density, ta, tb, wa, wb);
    end
    
    problems{t} = struct('theta', theta, 'W', W);
        
    mrf = fdReduction(theta, W);
    mrf2uai(mrf, sprintf('mrf_%d.UAI.LG', t));
    fprintf(1, 'Printed problem %d\n', t);
end
    
%     try
%       %% Our BBP
%       tic;
%       [logZ, oneMarg(:,t,1), twoMarg, misc] = BetheApprox_debug_mex(theta, W, epsilon, opts);        
%       betheTimes(t,1) = toc;    
%       bkEdges(t,1) = misc.nBKEdges;
%       maxFlowTimes(t,1)  = misc.BKMaxFlowTime;
%       intervalSzs(t,1) = misc.intervalSz;
%       
%       A(:,t,1) = misc.A;
%       B(:,t,1) = misc.B;
%       
%       fprintf(1, 'nBKNodes = %d; nBKEdges = %d\n', misc.nBKNodes, misc.nBKEdges);
%       times = [misc.makeMinSumTime misc.BKConstructTime misc.BKMaxFlowTime] ./ misc.mexTotTime;
%       
%       fprintf(1, 'Total time = %g; Fractions: Make minsum = %g, BK construction = %g, Max flow = %g, Rest = %g\n', ...
%               betheTimes(t), times(1), times(2), times(3), 1 - sum(times));
%       fprintf(1, 'MEX-call overhead fraction: %g\n', 1 - misc.mexTotTime / betheTimes(t));
%       
%       ABgap(:,t,1) = 1 - B(:,t,1) - A(:,t,1);
%       fprintf(1, '[A 1-B] gap mean = %g, max = %g, intervalSz = %g\n', mean(ABgap(:,t,1)), max(ABgap(:,t,1)), misc.intervalSz);
% 
%       %% Mooij bounds
%       tic;
%       [mooijlogZ, oneMarg(:,t,2), twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, mooijOpts);        
%       betheTimes(t,2) = toc;    
%       bkEdges(t,2) = misc.nBKEdges;
%       maxFlowTimes(t,2)  = misc.BKMaxFlowTime;
%       intervalSzs(t,2) = misc.intervalSz;
%       
%       A(:,t,2) = misc.A;
%       B(:,t,2) = misc.B;
%       
%       fprintf(1, 'nBKNodes = %d; nBKEdges = %d\n', misc.nBKNodes, misc.nBKEdges);
%       times = [misc.makeMinSumTime misc.BKConstructTime misc.BKMaxFlowTime] ./ misc.mexTotTime;
%       
%       fprintf(1, 'Total time = %g; Fractions: Make minsum = %g, BK construction = %g, Max flow = %g, Rest = %g\n', ...
%               betheTimes(t), times(1), times(2), times(3), 1 - sum(times));
%       fprintf(1, 'MEX-call overhead fraction: %g\n', 1 - misc.mexTotTime / betheTimes(t));
%       
%       ABgap(:,t,2) = 1 - B(:,t,2) - A(:,t,2);
%       fprintf(1, '[A 1-B] gap mean = %g, max = %g, intervalSz = %g\n', mean(ABgap(:,t,2)), max(ABgap(:,t,2)), misc.intervalSz);
% 
%       if runJT
%           % TODO: Consult other updates?
%           [trueLogZ, ~, trueOneMarg(:,t), trueTwoMarg, JTTimes(t)] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');           
%       end
%       
%       tic;
%       [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, lbpTimes(t)] = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');
%       
%       if checkMatlab
%           %% Check against MATLAB (cell execution)
%           [CHECKlogZ, CHECKoneMarg, CHECKtwoMarg, CHECKmisc] = solveBethe(theta, W, epsilon);
%           [logZ, oneMarg, twoMarg, misc] = BetheApprox_debug_mex(theta, W, epsilon, opts);
%           intervalSz = getIntervalSz(misc.A, misc.B, W, epsilon);
% 
%           assertElementsAlmostEqual(logZ, CHECKlogZ);    
%           assertElementsAlmostEqual(oneMarg, CHECKoneMarg);
%           assertElementsAlmostEqual(twoMarg, CHECKtwoMarg);
%           assertMiscEqual(misc, CHECKmisc);
%       end
%       
%       lbpGap  = lbpLogZ - logZ;
%       if lbpGap < 0
%           warning('LBP logZ - Approximate Bethe logZ gap was negative');
%       end
% 
%       lbpOneMargErr = mean(abs(lbpOneMarg - oneMarg(:,t)));
%       disp(['LBP logZ Gap = ' num2str(lbpGap)]);               
%       disp(['Average lbp-Bethe deviation of one-marginals: ' num2str(lbpOneMargErr) ]);
%       disp(['Average lbp-Bethe deviation of two-marginals: ' num2str(mean(abs(lbpTwoMarg(:) - twoMarg(:)))) ]);    
%       lbpLogZGaps(t) = lbpGap;
%       lbpTotDiff(t) = lbpOneMargErr;
%       
%       if abs(lbpGap) > epsilon
%           fn = ['questionable_case ' datestr(now, 0);];
%           warning(['LBP logZ - Approximate Bethe logZ gap exceeded epsilon. ' ...
%                    'Either LBP converged to a non-optimal stationary point ' ...
%                    'or the algorithm has a bug. Example saved as in' fn]);                 
%           save(fn);
%       end          
% 
%       disp(['Optimization finished with epsilon = ' num2str(epsilon)]);
%       
%       if runJT
%           trueGap = trueLogZ - logZ;
%           assert(trueGap >= 0, 'Approximate Bethe logZ did not lower bound true logZ');
%           assert(trueGap >= 0, 'Mooij Bethe logZ did not lower bound true logZ');
%           
%           oneMargErr = mean(abs(trueOneMarg(:,t) - oneMarg(:,t)));
%           
%           disp(['True logZ Gap = ' num2str(trueGap)]);
%           disp(['Average error of one-marginals: ' num2str(oneMargErr) ]);
%           disp(['Average error of two-marginals: ' num2str(mean(abs(trueTwoMarg(:) - twoMarg(:)))) ]);
%           trueLogZGaps(t)  = trueGap;
%           totDiff(t) = oneMargErr;
%       end
%     catch err
%       disp(err)
%     end
% 
%     save(fn); 
% end
% 
% betheTimes
% lbpTimes
% if runJT
%     JTTimes
% end
