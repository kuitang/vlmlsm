function [ theta, W ] = makeUnifProblem( nNodes, densityOrPattern, ta, tb, wb )
% makeProblem Make unbiased MRF with theta \in [ta, tb], W_ij \in [0, wb]
%
% [ theta, W ] = makeUnifProblem( nNodes, ta, tb, wa, wb )
%
% densityOrPattern - if scalar, then density; if sparse matrix, then
% pattern. Because W is symmetrized, this just lower-bounds the sparsity.
%
% theta - nNodes x 1
% W     - nNodes x nNodes graph with given sparsity, augmented to be
%         connected

    if issparse(densityOrPattern)
        W = sprand(densityOrPattern);
    else
        W = sprand(nNodes, nNodes, densityOrPattern);
    end    

    W = wb * W;                            % Scale and transform to [ta tb]
    W(1:nNodes+1:nNodes*nNodes) = 0;                      % Remove diagonal
    W = .5 * (W + W');                                         % Symmetrize                

    theta = unifrnd(ta, tb, nNodes, 1);
    theta = theta - 0.5 * sum(W, 2);                          % Remove bias
end

