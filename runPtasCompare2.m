%% Global (fixed) variables 

% Set dw so hard that the problems take nontrivial time.
dw = 6;

% Load the graphs etc.
fn       = 'PtasCompare_3.mat'
problems = cell(6, 100);

for nNodes = 4:10
    b = 2 * floor(log(nNodes));
    infn = sprintf('bregular%d.%d.mat', nNodes, b);
    load(infn);
    
    nGraphs = length(graphs);
    w = dw / b;

    for n = 1:10
        W     = w * sparse(graphs{n});
        theta = -0.5 * sum(W, 2);

        % Calling sequence
        fprintf(1, 'Starting bbp');

        tic;
        [bbpLogZ, ~, ~, bbpMisc] = BetheApprox_opt_mex(theta, W, 0.01, struct('useMooij', false));
        fprintf(1, 'nNodes = %d, dw = %d, bbp time = %g\n', nNodes, dw, toc);

        fprintf(1, 'Starting mooij');
        tic;
        [mkLogZ, ~, ~, mkMisc] = BetheApprox_opt_mex(theta, W, 0.01, struct('useMooij', true));
        fprintf(1, 'nNodes = %d, dw = %d, bbp time = %g\n', nNodes, dw, toc);

        fprintf(1, 'Starting brute');
        tic;
        [rawLogZ, ~, ~, rawMisc] = BetheApprox_opt_mex(theta, W, 0.01, struct('useMooij', false, 'maxIter', 0));
        fprintf(1, 'nNodes = %d, dw = %d, raw time = %g\n', nNodes, dw, toc);

        problems{nNodes, n} = var2struct(rawLogZ, bbpLogZ, mkLogZ, rawMisc, bbpMisc, mkMisc);
    end

    save(fn);
end

