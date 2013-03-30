function [ problems ] = unifLineSearch(j, params)
% unifLineSearch - Given an eta, search the J parameter space for bad LBP
%                  examples.
%
% eta      - a scalar
% params   - struct of nNodes, etaMin, etaMax, nPts; margThresh,
%            betheThresh, nSeqRnd
% problems - struct array of problem examples (oneMarg exceeds margThresh).
%            Contains fields
%             - eta, w
%             - lbpOneMarg, trueOneMarg, lbpTwoMarg, trueTwoMarg

    nNodes = params.nNodes;
    Wones  = ones(nNodes, nNodes);
    Wones  = Wones - diag(diag(Wones));
    W      = sparse(4*j * Wones);
    Nr = params.nSeqRnd;
    tones  = ones(nNodes, 1);    
    
    function tf = wrong(trueNegBethe, t1M, approxLogZ, a1M)
        % An instance is wrong iff the one marginals are far (absolute
        % difference > params.margThresh) and the Bethe free energies of the
        % approximate and true solutions are also far (absolute differences
        % > params.betheThresh).        
        tf = false;
        
        % Recall that for approximate solutions, Bethe = -logZ.        
        if trueNegBethe - approxLogZ > params.betheThresh && ...
            any(abs(t1M - a1M)) > params.margThresh
                tf = true;
                return;
        end        
    end    
    
    es = linspace(params.etaMin, params.etaMax, params.nPts);    
    nProblems = 0;
    
    for i = 1:params.nPts
        t = 2*es(i) - 2*(nNodes - 1)*j;
        theta = t * ones(nNodes, 1);                      
        
        % TODO: Figure out the initialization conditions in the codes.
        [trueLogZ, ~, trueOneMarg, trueTwoMarg, ~] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');
        [paLogZ, ~, paOneMarg, paTwoMarg, ~]       = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');        
        [sfLogZ, ~, sfOneMarg, sfTwoMarg, ~]       = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');                       
        
        trueNegBethe = -bethe(theta, W, trueOneMarg, trueTwoMarg);
                
        paFail = wrong(trueNegBethe, trueOneMarg, paLogZ, paOneMarg);
        sfFail = wrong(trueNegBethe, trueOneMarg, sfLogZ, sfOneMarg);        
        
        nWrong = 0;
        srBethes(params.nSeqRnd) = 0;
        srLogZs(params.nSeqRnd) = 0;
        
        % SeqRnd is harder; because it's random, we have to run many
        % trials.        
        srOneMargs(nNodes,Nr) = 0;
        for r = 1:Nr
            [srLogZs(r), ~, srOneMargs(:,r), ~, ~] = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');            
            if wrong(trueNegBethe, trueOneMarg, srLogZs(r), srOneMargs(:,r))
                nWrong = nWrong + 1;
            end
        end
        
        srFail = nWrong / Nr > 0.5;
        
        % Detect bifurcation (convergence to multiple fixed points).
        % However, due to random numerics, each convergence is slightly,
        % different, so we round. 6 digits after the decimal place seems
        % good enough, particularly since our tolerance was 1e-9 above.        
        [uniqSrLogZs, srIdxs] = unique(roundn(srLogZs, -6));
        % Filter down
        uniqSrOneMargs = srOneMargs(:,srIdxs);        
        
        % We must fail all three categories. Otherwise, some lower fixed
        % point does exist.
        if paFail && sfFail
            if ~srFail
                warning('SEQRND did not fail');
            end
            
            % Finally, test for feasibility.
            [A, B, alpha] = BBP(theta, W, 0.002, 10000);
            isz  = getIntervalSz(A, B, W, params.epsilon);
            nIntervals = sum((1 - B - A) ./ isz);
            if nIntervals > params.maxIntervals      
                fprintf(1, 'Infeasible problem at eta = %g, j = %g requires %d intervals.\n', ...
                        eta, j, nIntervals);
            else
                nProblems = nProblems + 1;            
                eta = es(i);
                problems(nProblems) = ...
                    var2struct(theta, W, eta, j, trueOneMarg, sfOneMarg, paOneMarg, uniqSrOneMargs, ...
                               trueLogZ, trueNegBethe, paLogZ, sfLogZ, uniqSrLogZs, ...
                               A, B, isz, nIntervals);
                problems(nProblems)
            end
        end
    end
    
    % Prevent error
    if nProblems == 0
        problems = false;
    end
end


