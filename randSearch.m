function [ problems, failVec ] = randSearch(params)
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
    
    %wrong = makePotentialWrong(params); 
    wrong = makeWrong(params); 
            
    nProblems = 0;
    failVec(3) = 0; 
    
    dstr = datestr(now, 0);
    fn = sprintf('randsearch_%d_nodes_%d_%s', params.nNodes, labindex, dstr); 
    fnFinal = sprintf('randsearch_%d_%s', labindex, dstr);
    
    % Print something and save intermediate results every 10%    
    decile = params.nIters / 10;
    countdown = decile;    

    for i = 1:params.nIters
        J   = params.jMax * sprand(params.adj);
        % Symmetrize!
        J   = 0.5 * (J + J'); 

        eta = unifrnd(params.etaMin, params.etaMax, params.nNodes, 1);
        
        [problem, fails] = searchCore(wrong, eta, J, params);
        failVec = failVec + fails;
        if ~isempty(fieldnames(problem))
            if ~params.quiet
                fprintf('Problem detected; totL1 = %g, nIntervals = %g\n', problem.totL1, problem.nIntervals);
            end

            nProblems = nProblems + 1;
            problems(nProblems) = problem;
        end
        
        % Periodically save and report progress.
        countdown = countdown - 1;
        if countdown == 0
            % Thin down the "problems" to the biggest ones before saving to disk
            if nProblems > 0
                [~, fastIdxs] = sort([problems.nIntervals]);
                fastProblems  = problems(fastIdxs(1:min(end, 5)));

                [~, idxs] = sort([problems.totL1], 1, 'descend');
                problems  = problems(idxs(1:min(end, 5)));
            end

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

    save(fnFinal);
end
