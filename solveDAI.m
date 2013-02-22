function [ logZ, joint, oneMarginals, twoMarginals ] = solveDAI(theta, W, method, daiOpts)
% solveDAI Find marginals of problem in Eq 1 by libDAI's methods
%   [oneMarginals, twoMarginals, joint, logZ] = solveJTree(theta, W)
%
%   theta - unary potentials
%   W     - (Sparse) pairwise interactions; only upper triangular taken
%
%   logZ         - exact (JTREE only) or approximate log partition function
%   joint        - final joint distribution
%   oneMarginals - nNodes x 1 vector of P(x_n = 1)
%   twoMarginals - 2 x 2 x nEdges array of 2x2 matrices where
%                  M(qi,qj) = P(x_i = qi, x_j = qj).
%                  We only compute pairwise marginals for edges that appear
%                  in W.

    [siVec, sjVec, swVec] = findUT(W);
    vars = [siVec sjVec];
    nNodes = length(theta);
    nEdges = length(siVec);
    
    % Convert the energy functional to a factor graph.
    psi     = cell(nNodes + nEdges, 1);
    varsets = cell(nNodes + nEdges, 1);
    for n = 1:nNodes        
        psi{n} = struct('Member', n, 'P', [1 exp(theta(n))]');
        varsets{n} = n;
    end
    
    for ne = 1:nEdges        
        w = swVec(ne);        
                
        P = ones(2,2);
        P(2,2) = exp(w);
                        
        psi{ne + nNodes} = struct('Member', vars(ne,:), 'P', P);
        varsets{ne + nNodes} = vars(ne,:);
    end
    
    % Call DAI
    [logZ,q,md,oneMStruct,twoMStruct,qmap] = dai(psi, method, daiOpts);

    % Sanitize the results
    joint = q{1}.P;
    % oneMStruct is ordered by the nodes because that's how we entered them.
    % At least, we hope so.
    oneMarginals = zeros(nNodes, 1);
    for n = 1:nNodes
        assert(oneMStruct{n}.Member == n);
        oneMarginals(n) = oneMStruct{n}.P(2);
    end
    
    twoMarginals = zeros(2, 2, nEdges);
    for ne = 1:nEdges
        mStructIdx = ne + nNodes;
        assert(all(twoMStruct{mStructIdx}.Member == vars(ne,:)));
        twoMarginals(:,:,ne) = twoMStruct{mStructIdx}.P;
    end
    
end

