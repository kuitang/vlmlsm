function [ psi ] = makePsi(theta, W)
% makePsi  Make potential structure (psi) for fastSolveDAI.
%   [ psi ] = makePsi(theta, W)
%
%   theta - unary potentials
%   W     - (Sparse) pairwise interactions; only upper triangular taken
%
%   psi   - structure to give to fastSolveDAI or dai

    [siVec, sjVec, swVec] = findUT(W);
    vars = [siVec sjVec];
    nNodes = length(theta);
    nEdges = length(siVec);
    
    % Convert the energy functional to a factor graph.
    psi     = cell(nNodes + nEdges, 1);    
    for n = 1:nNodes        
        psi{n} = struct('Member', n, 'P', [1 exp(theta(n))]');
        varsets{n} = n;
    end
    
    for ne = 1:nEdges        
        w = swVec(ne);        
                
        P = ones(2,2);
        P(2,2) = exp(w);
                        
        psi{ne + nNodes} = struct('Member', vars(ne,:), 'P', P);
    end    
end

