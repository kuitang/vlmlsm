% RUN ON A SEPARATE MACHINE FOR EACH NNODES
% set dw so that the problems take nontrivial time

dw = 6;

outfn = sprintf('VaryingRegular_%d_RESULTS.mat', nNodes);

infn = sprintf('varyingregular%d.mat', nNodes);
load(infn);

nGraphs = length(graphs);

results = cell(nGraphs, 1);

for j = 1:nGraphs
    G     = graphs{j};
    sG    = sum(graphs{j}, 2);
    d     = sG(1);
    w     = dw / d;

    W     = w * sparse(graphs{j});
    theta = -0.5 * sum(W ,2);

    tic;
    [mkLogZ, mkOneMarg, mkTwoMarg, mkMisc] = BetheApprox_opt_mex(theta, W, 0.01, struct('useMooij', true));
    mkTime = toc;
    fprintf(1, 'nNodes = %d, d = %d, dw = %g, mk time = %g\n', nNodes, d, dw, mkTime);

    results{j} = var2struct(nNodes, d, theta, W, mkLogZ, mkOneMarg, mkTwoMarg, mkMisc, mkTime);
    save(outfn);
end

