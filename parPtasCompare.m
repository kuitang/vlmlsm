function problems = parPtasCompare(graphs, w)
    parfor n = 1:12
        W     = w * sparse(graphs{n});
        nNodes = size(graphs{n}, 1);
        theta = -0.5 * sum(W, 2);

        % Calling sequence
        disp('Starting BBP');

        tic;
        [bbpLogZ, ~, ~, bbpMisc] = BetheApprox_opt_mex(theta, W, 0.01, struct('useMooij', false));
        fprintf(1, 'nNodes = %d, w = %g, bbp time = %g\n', nNodes, w, toc);

        disp('Starting MK');
        tic;
        [mkLogZ, ~, ~, mkMisc] = BetheApprox_opt_mex(theta, W, 0.01, struct('useMooij', true));
        fprintf(1, 'nNodes = %d, w = %g, mk time = %g\n', nNodes, w, toc);

        disp('Starting Brute');
        tic;
        [rawLogZ, ~, ~, rawMisc] = BetheApprox_opt_mex(theta, W, 0.01, struct('useMooij', false, 'maxIter', 0));
        fprintf(1, 'nNodes = %d, w = %d, brute time = %g\n', nNodes, w, toc);

        problems{n} = var2struct(theta, W, rawLogZ, bbpLogZ, mkLogZ, rawMisc, bbpMisc, mkMisc);
    end
end

