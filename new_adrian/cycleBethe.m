function [mu, lam] = cycleBethe(theta, W, varargin)
% [mu, lam] = cycleBethe(theta, W) minimize Bethe free energy over CYCLE.

    p = inputParser;
    p.addRequired('theta', @isnumeric);
    p.addRequired('W');
    p.addParamValue('epsilon', 0.1); % It's gonna be slow...        
    p.addParamValue('gapTol', 1e-5);
    p.addParamValue('initStep', 1e-2);
    p.parse(theta, W, varargin{:});
    
    o = p.Results;
    
    N = length(theta);
    W = triu(W, 1);
    assert(nnz(W) == N * (N - 1) / 2, 'W was not fully connected.');        
    
    % Chordless cycles in fully connected graphs are just triples. So
    % cycles is an Ncycles x 4 matrix, each row contains the nodes in the
    % cycle. (Columns 1 and 4 are repeated, to keep subsequent code simple)
            
    triplets  = combnk(1:N, 3);
    Ntriplets = size(triplets, 1);
    cycles(Ntriplets,6) = 0;
    for c = 1:length(triplets)
        edges = combnk(triplets(c,:), 2);
        cycles(c,:) = vec(edges')';
    end    
    
    % Make our constraint set
    Nconstr = 2 * N * (N - 1) * (N - 2) / 3;    
    
    % Each cycle has 4 Lagrange multipliers, due to the 4 different choices
    % for F.
    lam = ones(1, Nconstr);
    
    % Solve the initial problem without constraints.
    iter = 0;    
    
    % The dual should increase. Count the times it decreased instead. We
    % use the quantity to adjust stepsize [Rush12].    
    dualDecreaseIters = 0;    
    dual = 0;
    stepSz = o.initStep;
    oldDual = 0;
    viol = Inf;
    Econst = 0;
    
    % DATA STRUCTURE TO STORE BOTH CHORDLESS CYCLES AND THEIR ODD EDGE
    % SUBSETS?
    %
    % Odd edge subsets have variable size.
    
    % Helper for the loop
    function reparam(edge, dTheta, dW)
        theta(edge) = theta(edge) + dTheta;

        u = edge(1); v = edge(2);
        W(u,v) = W(u,v) + dW;
    end

    while abs(viol) > o.gapTol                
        % Solve \argmin_{\mu} L(\mu, \lambda).
        %
        % Absorb constraint terms into the primal
        for n = 1:Ntriplets           % for all chordless cycles
            for i = 1:3               % for all odd edge subsets of size 1
                cyc      = cycles(n,:);
                oddEdge  = cyc(i:(i+1));
                compEdges = cyc;
                compEdges(i:(i+1)) = [];
                
                % Each cycle has 4 constraints.
                iconstr  = 4*(n-1) + i;
                lm       = lam(iconstr);
                
                % Reparameterize the odd edge
                reparam(oddEdge, lm, -2*lm);                
                
                % Econst has to update because we updated parameters
                % asymmetrically; see Algorithm 1 in Adrian's notes.
                Econst = Econst - lm;
                
                % Reparameterize the complement edges (there are two)
                for k = 1:2:4
                    ce = compEdges(k:(k+1));                    
                    reparam(ce, -lm, 2*lm);
                    
                    % No change the Econst                    
                end
            end
            
            % Reparameterize the chordless cycle with all 3 edges in F
            iconstr = 4*n;
            lm = lam(iconstr);
            for k = 1:2:6
                oddEdge = cyc(k:(k+1));
                reparam(oddEdge, lm, -2*lm);
                
                Econst = Econst - lm;
            end
        end

        % Solve it.
        [mu, xi, primal] = solveBetheExact(theta, W, 'epsilon', o.epsilon);
        
        % We want a symmetric xi
        xi = triu(xi, 1);
        xi = xi + xi';
        
        % Solve \argmax_{\lambda} L(\mu, \lambda) by subgradient ascent.
        
        % Adrian's decomposition of the term in the parentheses.
        %
        % The subgradient is just the stuff inside the constraint
        % parentheses--the violations.
        muMat    = bsxfun(@plus, mu, mu');
        mainDiag = 1 + 2*xi - muMat;
        offDiag  = muMat - 2*xi;
        
        subGrad(Nconstr) = 0;
        
        constr = 1;
        for n = 1:Ntriplets           % for all chordless cycles
            for i = 1:3               % for all odd edge subsets F of size 1
                cyc      = cycles(n,:);
                oddEdge  = cyc(i:(i+1));
                compEdges = cyc;
                
                % Summation calculations are quite inelegant.
                compEdges(i:(i+1)) = [];                
                
                uu = oddEdge(1); vv = oddEdge(2);
                is = compEdges([1,3]);
                js = compEdges([2,4]);                                
                compInds = sub2ind([N, N], is, js);
                
                subGrad(constr) = 1 - mainDiag(uu,vv) - sum(offDiag(compInds));
                
                constr = constr + 1;
            end
            
            % the chordless cycle with all 3 edges in F
            oddEdges = cycles(n,:);
            is = oddEdges(1:2:6);
            js = oddEdges(2:2:6);
            oddInds = sub2ind([N, N], is, js);
            
            subGrad(constr) = 1 - sum(mainDiag(oddInds));
            constr = constr + 1;            
        end        
        
        % Calculate statistics for this iteration, before stepping.
        oldDual = dual;
        viol = sum(lam .* subGrad);
        dual = primal + viol;
        
        assert(viol >= 0, 'the dual must lower-bound the primal!');
        
        % Step size schedule in in [Rush12]; see [Komodakis11] for other ideas.
        if oldDual > dual
            % We overstepped, so scale down our stepsize.
            dualDecreaseIters = dualDecreaseIters + 1;
            stepSz = o.initStep / (1 + dualDecreaseIters);
        end
        
        fprintf(1, 'Iter %d, primal = %g, dual violations (gap) = %g\n', iter, primal, viol);
               
        % lam is set for next iteration.
        
        % Projected subgradient ascent (multipliers must be positive), see
        % [Komodakis11].
        lam = max(0, lam + stepSz * subGrad);                
        iter = iter + 1;
        
        % Save the state somewhere???
    end
    
    fprintf(1, 'Duality gap < %g', o.gapTol);
end

