function randSparse(dw, T, graphs)
    runBBP = false;
    runMooij = true;

    nInnerTrials = 10;
    nGraphs = length(graphs);
    nTrials = nInnerTrials * nGraphs;

    %% Setup
    path_to_dai = '../libDAI-0.3.1/matlab';
    addpath(path_to_dai);

    %% Problem
    % For nNodes = 50, this runs out of memory :(

    nNodes = size(graphs{1}, 1);

    runJT = false;
    %testTrees = true;    

    epsilon = 0.01;

    oneMarg(nNodes,nTrials,2) = 0;
    trueOneMarg(nNodes,nTrials) = 0;
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

    makeMinSumTimes(nTrials,2) = 0;

    lbpTimes(nTrials) = 0;
    JTTimes(nTrials) = 0;

    opts = struct('useMooij', true);
    mooijOpts = struct('useMooij', true);

    problems = cell(1, nTrials);
        
    fn = sprintf('mk_bregular_rand_dw_%d_T_%d_nNodes_%d', dw, T, nNodes);        

    %% Loop it
    deg = sum(graphs{1}(:,1), 1);
    w   = dw / deg;
    
    for n = 1:nGraphs
        for it = 1:nInnerTrials
            t = nInnerTrials*(n - 1) + it;

            %% Set up the problem
            tt    = unifrnd(-T, T, nNodes, 1);
            % Important: symmetrize!
            W     = w * sprand(graphs{n});
            W     = 0.5 * (W + W');
            %W      = sprand(nNodes, nNodes, 0.1);
            %W      = 0.5 * (W + W');

            theta = tt - sum(W, 2);
            problems{t} = struct('theta', theta, 'W', W);

            %%
            try
                %% LBP        
                [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, lbpTimes(t)] = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');

                %% Other runs
                if runJT
                    % TODO: Consult other updates?
                    [trueLogZ, ~, trueOneMarg(:,t), trueTwoMarg, JTTimes(t)] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');
                end

                %% Our BBP
                if runBBP
                    tic;


                    [logZ, oneMarg(:,t,1), twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, opts);



                    betheTimes(t,1) = toc;
                    makeMinSumTimes(t,1) = misc.makeMinSumTime;
                    bkEdges(t,1) = misc.nBKEdges;
                    maxFlowTimes(t,1)  = misc.BKMaxFlowTime;
                    intervalSzs(t,1) = misc.intervalSz;

                    A(:,t,1) = misc.A;
                    B(:,t,1) = misc.B;

                    fprintf(1, 'nBKNodes = %d; nBKEdges = %d\n', misc.nBKNodes, misc.nBKEdges);
                    times = [misc.makeMinSumTime misc.BKConstructTime misc.BKMaxFlowTime] ./ misc.mexTotTime;

                    fprintf(1, 'Total time = %g; Fractions: Make minsum = %g, BK construction = %g, Max flow = %g, Rest = %g\n', ...
                        betheTimes(t), times(1), times(2), times(3), 1 - sum(times));
                    fprintf(1, 'MEX-call overhead fraction: %g\n', 1 - misc.mexTotTime / betheTimes(t));

                    ABgap(:,t,1) = 1 - B(:,t,1) - A(:,t,1);
                    fprintf(1, '[A 1-B] gap mean = %g, max = %g, intervalSz = %g\n', mean(ABgap(:,t,1)), max(ABgap(:,t,1)), misc.intervalSz);

                    %% LBP Errors        
                    lbpGap  = lbpLogZ - logZ;
                    if lbpGap < 0
                        warning('LBP logZ - Approximate Bethe logZ gap was negative');
                    end

                    lbpOneMargErr = mean(abs(lbpOneMarg - oneMarg(:,t,1)));
                    disp(['Bethe LBP logZ Gap = ' num2str(lbpGap)]);        

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
                end

                if runMooij
                    %% Mooij bounds
                    tic;



                    [mooijLogZ, oneMarg(:,t,2), twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, mooijOpts);




                    betheTimes(t,2) = toc;
                    makeMinSumTimes(t,2) = misc.makeMinSumTime;
                    bkEdges(t,2) = misc.nBKEdges;
                    maxFlowTimes(t,2)  = misc.BKMaxFlowTime;
                    intervalSzs(t,2) = misc.intervalSz;

                    A(:,t,2) = misc.A;
                    B(:,t,2) = misc.B;

                    fprintf(1, 'nBKNodes = %d; nBKEdges = %d\n', misc.nBKNodes, misc.nBKEdges);
                    times = [misc.makeMinSumTime misc.BKConstructTime misc.BKMaxFlowTime] ./ misc.mexTotTime;

                    fprintf(1, 'Total time = %g; Fractions: Make minsum = %g, BK construction = %g, Max flow = %g, Rest = %g\n', ...
                        betheTimes(t), times(1), times(2), times(3), 1 - sum(times));
                    fprintf(1, 'MEX-call overhead fraction: %g\n', 1 - misc.mexTotTime / betheTimes(t));

                    ABgap(:,t,2) = 1 - B(:,t,2) - A(:,t,2);
                    fprintf(1, '[A 1-B] gap mean = %g, max = %g, intervalSz = %g\n', mean(ABgap(:,t,2)), max(ABgap(:,t,2)), misc.intervalSz);

                    %% Mooij errors
                    mooijLbpGap  = lbpLogZ - mooijLogZ;

                    mooijLbpOneMargErr = mean(abs(lbpOneMarg - oneMarg(:,t,2)));
                    disp(['Mooij LBP logZ Gap = ' num2str(mooijLbpGap)]);        

                    mooijLbpLogZGaps(t) = mooijLbpGap; 
                    mooijLbpTotDiff(t)  = mooijLbpOneMargErr;

                    disp(['Average lbp-Bethe deviation of one-marginals: ' num2str(mooijLbpOneMargErr) ]);
                    disp(['Average lbp-Bethe deviation of two-marginals: ' num2str(mean(abs(lbpTwoMarg(:) - twoMarg(:)))) ]);


                end

                disp(['Optimization finished with epsilon = ' num2str(epsilon)]);

                if runJT
                    trueGap = trueLogZ - logZ;            

                    assert(trueGap >= 0, 'Approximate Bethe logZ did not lower bound true logZ');
                    assert(trueGap >= 0, 'Mooij Bethe logZ did not lower bound true logZ');

                    oneMargErr = mean(abs(trueOneMarg(:,t) - oneMarg(:,t)));

                    disp(['True logZ Gap = ' num2str(trueGap)]);
                    disp(['Average error of one-marginals: ' num2str(oneMargErr) ]);
                    disp(['Average error of two-marginals: ' num2str(mean(abs(trueTwoMarg(:) - twoMarg(:)))) ]);
                    trueLogZGaps(t)  = trueGap;
                    totDiff(t) = oneMargErr;
                end

                save(fn);

            catch err
                disp(err);
                problems{t} = err;
            end
        end

    end

    betheTimes
    lbpTimes
    if runJT
        JTTimes
    end

    save(fn)
end
