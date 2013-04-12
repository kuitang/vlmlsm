function [ problem, failVec ] = searchCore(wrong, eta, J, params)
% searchCore Evaluate an example (eta, J) for "wrongness".
%
    problem = struct();
    nNodes = size(J, 1);
    nEdges = nnz(J) / 2;    
    Nr     = params.nSeqRnd;
    
    % Convert to our parameterization
    theta = 2*eta - 2*sum(J, 2);
    W     = 4*J;

    psi = makePsi(theta, W);

    %[ JBack, etaBack ] = mooijParam(W, theta);
    %assertElementsAlmostEqual(JBack, J);
    %assertElementsAlmostEqual(etaBack, eta);

    % TODO: Figure out the initialization conditions in the codes.    
    [trueLogZ, trueOneMarg, ~, trueTwoMarg] = fastSolveDAI(nNodes, nEdges, psi, 'JTREE', '[updates=HUGIN,verbose=0]');
%         [checkLogZ, ~, checkOneMarg, checkTwoMarg, ~] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');
%         assertElementsAlmostEqual(checkLogZ, trueLogZ);
%         assertElementsAlmostEqual(trueOneMarg, checkOneMarg);
%         assertElementsAlmostEqual(trueTwoMarg, checkTwoMarg);
    [paLogZ, paOneMarg, paLogZHist] = fastSolveDAI(nNodes, nEdges, psi, 'BP', '[inference=SUMPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');        
    [sfLogZ, sfOneMarg, sfLogZHist] = fastSolveDAI(nNodes, nEdges, psi, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');                       

    if params.plot
        figure(1);
        plot(paLogZHist);
        title('PARALL logZ');

        figure(2);
        plot(sfLogZHist);
        title('SEQFIX logZ');
    end
    
    trueNegBethe = -bethe(theta, W, trueOneMarg, trueTwoMarg);

    paFail = wrong(trueNegBethe, trueOneMarg, paLogZ, paOneMarg);
    sfFail = wrong(trueNegBethe, trueOneMarg, sfLogZ, sfOneMarg);        

    paL1   = sum(abs(trueOneMarg - paOneMarg));
    sfL1   = sum(abs(trueOneMarg - sfOneMarg));

    % SeqRnd is harder; because it's random, we have to run many
    % trials.          
    nWrong = 0;        
    srLogZs(Nr) = 0;
    srOneMargs(nNodes,Nr) = 0;
    for r = 1:Nr
        [srLogZs(r), srOneMargs(:,r), ~] = fastSolveDAI(nNodes, nEdges, psi, 'BP', '[inference=SUMPROD,updates=SEQRND,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');            
        if wrong(trueNegBethe, trueOneMarg, srLogZs(r), srOneMargs(:,r))
            nWrong = nWrong + 1;
        end
    end
    srFail = nWrong / Nr > 0.5;    
    
    failVec = [paFail sfFail srFail];

    % Detect bifurcation (convergence to multiple fixed points).
    % However, due to random numerics, each convergence is slightly,
    % different, so we round. 6 digits after the decimal place seems
    % good enough, particularly since our tolerance was 1e-9 above.

    roundSrLogZs = roundn(srLogZs, -6);
    [uniqSrLogZs, srIdxs] = unique(roundSrLogZs);
    nFixedPoints = length(roundSrLogZs);
    distSrLogZs = histc(roundSrLogZs, uniqSrLogZs) ./ nFixedPoints;
    
    % Filter down
    uniqSrOneMargs = srOneMargs(:,srIdxs);        

    % We must fail all three categories. Otherwise, some lower fixed
    % point does exist.
    %
    % OR-FAIL (yadiyadiyah)
    if paFail || sfFail
        if ~srFail
            warning('SEQRND did not fail');
        end
        % Finally, test for feasibility.
        [A, B, ~] = BBP(theta, W, 0.0002, 10000);
        isz  = getIntervalSz(A, B, W, params.epsilon);
        nIntervals = sum((1 - B - A) ./ isz);        
        if nIntervals < 0
            warning('Got negative nIntervals; values probably too large');
        end

        if abs(nIntervals) > params.maxIntervals      
            %fprintf(1, 'Infeasible problem required %d intervals.\n', ...
            %        nIntervals);
        else
            totL1 = paL1 + sfL1;
            problem = ...
                var2struct(theta, W, eta, J, trueOneMarg, sfOneMarg, paOneMarg, uniqSrOneMargs, ...
                           trueLogZ, trueNegBethe, paLogZ, sfLogZ, uniqSrLogZs, distSrLogZs, ...
                           A, B, isz, nIntervals, paL1, sfL1, totL1);            
        end
    end
end

