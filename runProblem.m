function runProblem( problem, epsilon )
% runProblem - Run LBP, JTree, and our algorithm for one identified problem

    theta = problem.theta;
    W     = problem.W;

    opts.useMooij = false;    
    [bethelogZ, betheOneMarg, betheTwoMarg, misc]    = BetheApprox_debug_opt(theta, W, epsilon, opts);
    [trueLogZ, ~, trueOneMarg, trueTwoMarg, JTTimes] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');           
    [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, lbpTimes]   = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');

end

