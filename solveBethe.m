function [logZ, oneMarginals, twoMarginals, misc] = solveBethe(theta, W, epsilon)
% solveBethe Solve a binary MRF problem with the Bethe bound approximation
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
    
    [misc.A, misc.B, alpha] = BBP(theta, W);
    % separate out for easy debugging
    misc.intervalSz = getIntervalSz(misc.A, misc.B, W, epsilon);
    [misc.D, WW, Vi, misc.Vm, qr] = boundMRF(theta, W, misc.A, misc.B, misc.intervalSz);
    [x, e, misc.elMat] = MultiLabelSubModular(misc.D, WW, Vi, misc.Vm);
    misc.e = e;
    
    % How big was our problem?
    elNzRows = (misc.elMat(:,1) ~= 0) & (misc.elMat(:,2) ~= 0);    
    misc.nDiscreteEdges = sum(elNzRows);
    
    % e is a vector of energies; e(1) is total energy.
    logZ = -e(1);

    % Translate levels back to marginals
    oneMarginals = zeros(nNodes, 1);
    for n = 1:nNodes
        oneMarginals(n) = qr{n}(x(n));
    end
  
    twoMarginals = zeros(2, 2, nEdges);
    for ne = 1:nEdges
        i = siVec(ne);
        j = sjVec(ne);

        a = alpha(i,j);
        q_i = oneMarginals(i);
        q_j = oneMarginals(j);

        twoMarginals(:,:,ne) = marginalize(a, q_i, q_j);
    end
    
    for ne = 1:nEdges
        i = siVec(ne);
        j = sjVec(ne);

        checkIMarginal = sum(twoMarginals(:,:,ne), 2);   
        checkJMarginal = sum(twoMarginals(:,:,ne), 1);

        assertElementsAlmostEqual(oneMarginals(i), checkIMarginal(2));
        assertElementsAlmostEqual(oneMarginals(j), checkJMarginal(2));
    end

end

