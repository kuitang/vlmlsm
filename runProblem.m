function runProblem( problemPath, userEpsilon )
% runProblem - Run LBP, JTree, and our algorithm for one identified problem

    path_to_dai = '../libDAI-0.3.1/matlab';
    addpath(path_to_dai);

    load(problemPath);
    problem = problems(1);

    % Maybe?
    epsilon = userEpsilon;
    problem.epsilon = userEpsilon;

    theta = problem.theta;
    W     = problem.W;

    opts.useMooij = false;    
    [bethelogZ, betheOneMarg, betheTwoMarg, misc]    = BetheApprox_opt_mex(theta, W, userEpsilon, opts);
    [trueLogZ, ~, trueOneMarg, trueTwoMarg, JTTimes] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');           
    [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, lbpTimes]   = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');

    disp(['Optimization finished with epsilon = ' num2str(userEpsilon)]);

    savePath = [problemPath '_RESULTS_NEW'];
    save(savePath);

end

