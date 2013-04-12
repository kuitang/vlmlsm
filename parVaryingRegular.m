function results = parVaryingRegular(graphs, dw)
% We have to do it this way to avoid a transparency error

    nGraphs = length(graphs);
    results = cell(nGraphs, 1);

    parfor j = 1:nGraphs
        G      = graphs{j}
        nNodes = size(G, 1);
        sG     = sum(G, 2);
        d      = sG(1);
        w      = dw / d;

        W     = w * sparse(G);
        theta = -0.5 * sum(W ,2);

        tic;
        [mkLogZ, mkOneMarg, mkTwoMarg, mkMisc] = BetheApprox_opt_mex(theta, W, 0.01, struct('useMooij', true));
        mkTime = toc;
        fprintf(1, 'nNodes = %d, d = %d, dw = %g, mk time = %g\n', nNodes, d, dw, mkTime);

        results{j} = var2struct(nNodes, d, theta, W, mkLogZ, mkOneMarg, mkTwoMarg, mkMisc, mkTime);
    end

end

