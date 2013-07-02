function [logZ, oneMarginals, twoMarginals, misc] = solveBetheNew(theta, W, gams)
% solveBethe Solve a binary MRF problem with given gams
%
%   [logZ, oneMarginals, twoMarginals] = solveBethe(theta, W);
%   theta   - 1 x nNodes unary potentials
%   W       - (Sparse) pairwise interactions; only upper triangular kept.
%             nEdges nonzero entries in the upper triangular. (Should be
%             submodular; results not guaranteed otherwise.)
%   epsilon - bound on difference to true logZ (negative free energy)
%   
%   logZ  - Log partition function, which is equal to the negative minimum
%           Bethe free energy.
%   oneMarginals - 1 x nNodes vector of P(X_n = 1)
%   twoMarginas  - 2 x 2 x nEdges matrix of P(X_i, X_j).
%   misc - structure of diagnostic information

    nNodes = length(theta);
    [siVec, sjVec, swVec] = findUT(W);    
    nEdges = length(siVec);        
    assert(size(W, 1) == nNodes);
    
    % Symmetrize W (BBP and boundMRP expect symmetric W)
    W = sparse(vertcat(siVec, sjVec), ...
               vertcat(sjVec, siVec), ...
               vertcat(swVec, swVec), ...
               nNodes, nNodes);
    
    assert(all(W(:) >= 0), 'W contains negative entries. Not submodular!');    
           
    [D, newW, Vi, Vm] = boundMRFNew(theta, W, gams);
    [x, e]  = MultiLabelSubModular(D, newW, Vi, Vm);
    
    misc = var2struct(D, newW, Vi, Vm, x, e);
    
    logZ = -e(1);

    % Translate levels back to marginals
    oneMarginals = zeros(nNodes, 1);
    for n = 1:nNodes
        oneMarginals(n) = gams{n}(misc.x(n));
    end
  
    twoMarginals = zeros(2, 2, nEdges);
    alpha = exp(abs(W)) - 1;
    for ne = 1:nEdges
        i = siVec(ne);
        j = sjVec(ne);

        a = alpha(i,j);
        q_i = oneMarginals(i);
        q_j = oneMarginals(j);

        twoMarginals(:,:,ne) = reshape(marginalize(a, q_i, q_j), 2, 2);
    end
    
    for ne = 1:nEdges
        i = siVec(ne);
        j = sjVec(ne);

        checkIMarginal = sum(twoMarginals(:,:,ne), 2);   
        checkJMarginal = sum(twoMarginals(:,:,ne), 1);

        assertElementsAlmostEqual(oneMarginals(i), checkIMarginal(2), 'relative', 1e-5);
        assertElementsAlmostEqual(oneMarginals(j), checkJMarginal(2), 'relative', 1e-5);
    end

end

