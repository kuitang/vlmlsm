function [intervalSz, varargout] = getIntervalSz(A, B, W, epsilon)
% getIntervalSz Calculate the required interval size for epsilon accuracy
%
% intervalSz = getIntervalSz(A, B, W, epsilon)    
% [intervalSz, kappa, theoryBound, aMax, bMax] = getIntervalSz(A, B, W, epsilon)
    nNodes = size(W, 1);
    nEdges = nnz(W) / 2;
    
    alpha = exp(abs(W)) - 1;
    [iVec, jVec, aVec] = findUT(alpha);
    
    deg = zeros(nNodes, 1);
    for n = 1:nNodes
        deg(n) = sum(W(:,n) > 0);
        assert(deg(n) > 0, 'All nodes must be connected');
    end
    
    % (Eq 16)
    eta = min(A, B);
    % (Theorem 15)
    aMax = -Inf;
    bMax = -Inf;
    for ne = 1:nEdges
        i = iVec(ne);
        j = jVec(ne);
        
        aij = aVec(ne);
        ei = eta(i);
        ej = eta(j);
        
        aa = aij*(aij + 1) / (4 * (2*aij + 1) * ei * ej * (1 - ei) * (1 - ej));
        if aa > aMax
            aMax = aa;
        end        
    end
    
    for i = 1:nNodes
        neighborAlphas = alpha(i,W(:,i) > 0);
        
        bInner = 1 - deg(i) + ...
                 sum( (neighborAlphas + 1).^2 ./ ...
                      (2*neighborAlphas + 1) );                
        bb = 1 / (eta(i) * (1 - eta(i))) * bInner;
        if bb > bMax
            bMax = bb;
        end
    end
    
    % (Eq 18)
    Omega = max(aMax, bMax);

    % Use the stupider bound    
    density = (nnz(W) + nNodes) / numel(W);           % \Sigma in the paper
    Lambda = nNodes * Omega * sqrt(density);                                   
    if nargout >= 2              % output kappa (\Sigma^{3/4} \Omega^{3/2})
        varargout{1} = density^(3/4) * Omega^(3/2);
    end
    if nargout >= 3 % output kappaBound 
        varargout{2} = (nNodes^6 * Omega^(3/2) * density^(3/4) / (sqrt(8) * epsilon^(3/2)));
    end
    if nargout >= 5
        varargout{3} = aMax;
        varargout{4} = bMax;
    end
    intervalSz = sqrt(2*epsilon / (Lambda));
    
end
