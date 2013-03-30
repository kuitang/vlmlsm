function [ problems, nFails ] = randSearch(params)
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
    
    nNodes = params.nNodes;
    nEdges = nnz(params.adj) / 2;
    Nr = params.nSeqRnd;
    nProblems = 0;
    
    nFails = struct('sf', 0, 'pa', 0, 'sr', 0);
    
    fn = ['randsearch_checkpoint' labindex '_' datestr(now, 0)];
    
    % Print something and save intermediate results every 10%    
    decile = params.nIters / 10;
    countdown = decile;    

    for i = 1:params.nIters
        J   = params.jMax * sprand(params.adj);
        eta = unifrnd(params.etaMin, params.etaMax, params.nNodes, 1);
        
        % Convert to our parameterization
        theta = 2*eta - 2*sum(J, 2);
        W     = 4*J;
        
        psi = makePsi(theta, W);
        
        %[ JBack, etaBack ] = mooijParam(W, theta);
        %assertElementsAlmostEqual(JBack, J);
        %assertElementsAlmostEqual(etaBack, eta);
        
        % TODO: Figure out the initialization conditions in the codes.
        [trueLogZ, trueOneMarg, trueTwoMarg] = fastSolveDAI(nNodes, nEdges, psi, 'JTREE', '[updates=HUGIN,verbose=0]');
%         [checkLogZ, ~, checkOneMarg, checkTwoMarg, ~] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');
%         assertElementsAlmostEqual(checkLogZ, trueLogZ);
%         assertElementsAlmostEqual(trueOneMarg, checkOneMarg);
%         assertElementsAlmostEqual(trueTwoMarg, checkTwoMarg);
        [paLogZ, paOneMarg]       = fastSolveDAI(nNodes, nEdges, psi, 'BP', '[inference=SUMPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');        
        [sfLogZ, sfOneMarg]       = fastSolveDAI(nNodes, nEdges, psi, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');                       
        
        trueNegBethe = -bethe(theta, W, trueOneMarg, trueTwoMarg);
                
        paFail = wrong(trueNegBethe, trueOneMarg, paLogZ, paOneMarg);
        if paFail
            nFails.pa = nFails.pa + 1;
        end
        
        sfFail = wrong(trueNegBethe, trueOneMarg, sfLogZ, sfOneMarg);        
        if sfFail
            nFails.sf = nFails.sf + 1;
        end
        
        % SeqRnd is harder; because it's random, we have to run many
        % trials.          
        nWrong = 0;        
        srLogZs(Nr) = 0;
        srOneMargs(nNodes,Nr) = 0;
        for r = 1:Nr
            [srLogZs(r), srOneMargs(:,r)] = fastSolveDAI(nNodes, nEdges, psi, 'BP', '[inference=SUMPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');            
            if wrong(trueNegBethe, trueOneMarg, srLogZs(r), srOneMargs(:,r))
                nWrong = nWrong + 1;
            end
        end
        
        srFail = nWrong / Nr > 0.5;
        if srFail
            nFails.sr = nFails.sr + 1;
        end
        
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
            eta = es(i);            
            if nIntervals > params.maxIntervals      
                fprintf(1, 'Infeasible problem at requires %d intervals.\n', ...
                        eta, j, nIntervals);
            else
                nProblems = nProblems + 1;                            
                problems(nProblems) = ...
                    var2struct(theta, W, eta, j, trueOneMarg, sfOneMarg, paOneMarg, uniqSrOneMargs, ...
                               trueLogZ, trueNegBethe, paLogZ, sfLogZ, uniqSrLogZs, ...
                               A, B, isz, nIntervals);
                problems(nProblems)
            end
        end
        
        % Periodically save and report progress.
        countdown = countdown - 1;
        if countdown == 0
           fprintf('%s -- Lab %d: Iter %d of %d: %d problems\n', ...
                   datestr(now, 0), labindex, i, params.nIters, nProblems);           
           save(fn);
           countdown = decile;
        end
    end
    
    % Prevent error: assign empty struct.
    if nProblems == 0
        problems = struct();
    end
end
