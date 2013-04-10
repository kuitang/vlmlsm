function worstSparse(dW, T, nNodes, epsilon)
    path_to_dai = '../libDAI-0.3.1/matlab';
    addpath(path_to_dai);

    epsilon = 0.1;

    opts = struct('useMooij', false);

    adj = diag(ones(nNodes - 1, 1), 1);
    adj(nNodes,1) = 1;

    % Now add extra edges
    deg = log(nNodes) / log(2);
    extraEdges = deg - 2;

    for n = 0:(nNodes-1)
        extraNodes = mod(n + 2*(1:extraEdges) + 1, nNodes) + 1;
        %[n + 1 extraNodes]
        adj(n+1,extraNodes) = 1;
    end

    % Symmetrize
    adj = adj | adj';

    W = sparse(dW / deg * adj)

    Tee = T * ones(nNodes, 1);
    Tee(2:2:nNodes) = -T;

    theta = (Tee - dW/2) .* ones(nNodes, 1)

    % Interval size calculation
    [A, B, alpha] = BBP(theta, W);
    [iSz, ~]      = getIntervalSz(A, B, W, epsilon);
    nIntervals    = sum(1 - B - A) / iSz

    [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, lbpTime] = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');
    [trueLogZ, ~, trueOneMarg, trueTwoMarg, JTTimes] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');
    [logZ, oneMarg, twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, opts);

    [misc.A 1 - misc.B]

    fn = sprintf('worst_dW_%d_T_%d_nNodes_%d\n', dW, T, nNodes);
    save(fn)
end

