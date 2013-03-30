function [ problems, failVec ] = unifLineSearch(j, params)
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
    
    wrong = makeWrong(params);

    J = sparse(j * params.adj);
    failVec(3) = 0;
    es = linspace(params.etaMin, params.etaMax, params.nPts);    
    nProblems = 0;
    
    for i = 1:params.nPts
        eta = es(i);
                
        [problem, fails] = searchCore(wrong, eta, J, params);
        failVec = failVec + fails;
        if ~isempty(fieldnames(problem))
            nProblems = nProblems + 1;
            problems(end+1) = problem;
        end 
    end

    % Prevent error
    if nProblems == 0
        problems = false;
    end
end


