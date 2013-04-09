function runProblem( problemPath, epsilon )
% runProblem - Run LBP, JTree, and our algorithm for one identified problem

    path_to_dai = '../libDAI-0.3.1/matlab';
    addpath(path_to_dai);

    load(problemPath);
    problem = problems(1);

    theta = problem.theta;
    W     = problem.W;

    opts.useMooij = false;    
    [bethelogZ, betheOneMarg, betheTwoMarg, misc]    = BetheApprox_opt_mex(theta, W, epsilon, opts);
    [trueLogZ, ~, trueOneMarg, trueTwoMarg, JTTimes] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');           
    [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, lbpTimes]   = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');

    disp(['Optimization finished with epsilon = ' num2str(epsilon)]);

    savePath = [problemPath '_RESULTS'];
    save(savePath);

end

