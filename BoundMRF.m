function [ D, newW, Vi, Vm ] = boundMRF( theta, W, A, B, alpha, epsilon )
% Construct an MRF for the Bethe bound approximation.
%
% Output an MRF with N nodes and, each having as many labels as
% discretization intervals determined by A(n) : (1 - B(n) : gamma.
%
% theta - N-vector of unary weights
% W     - N x N symmetric sparse matrix of pairwise edge weights (Eq 1)
% A     - N-vector of lower bounds
% B     - N-vector of complementary upper bounds
% alpha - output from BBP
% epsilon - tolerance from the true optimum value
%
% D, W, Vi, Vm - Arguments to pass to MultiLabelSubModular
%
% Equations and numbers taken from the 31 Dec 2012 draft.

    % From the symmetric graph, extract an edge list where iVec(n) < jVec(n)
    nNodes = length(theta);
    [iVec, jVec, wVec] = find(W);
    sel = iVec < jVec;
    iVec = iVec(sel);
    jVec = jVec(sel);
    wVec = wVec(sel);
    
    nEdges = length(iVec);    

    deg = zeros(nNodes);
    for n = 1:nNodes
        deg(n) = sum(W(:,n) > 0);
        assert(deg(n) > 0, 'All nodes must be connected');
    end
    
    % Interval size (gamma) calculation
    
    % (Eq 16)
    eta = min(A, B);
    % (Theorem 15)
    aMax = -Inf;
    bMax = -Inf;
    for ne = 1:nEdges
        i = iVec(ne);
        j = iVec(ne);
        
        aij = alpha(i,j);
        ei = eta(i);
        ej = eta(j);
        
        aa = aij*(aij + 1) / (4 * (2*aij + 1) * ei * ej * (1 - ei) * (1 - ej));
        if aa > aMax
            aMax = aa;
        end        
    end
    
    for i = 1:nNodes
        neighborAlphas = alpha(i,W(i,:) > 0);
        
        bInner = 1 - deg(ne) + ...
                 sum( (neighborAlphas + 1).^2 ./ ...
                      (2*neighborAlphas + 1) );                
        bb = 1 / (eta(i) * (1 - eta(i))) * bInner;
        if bb > bMax
            bMax = bb;
        end
    end
    
    % (Eq 18)
    Omega = max(aMax, bMax);
    density = nnz(W) / numel(W);
    eigenBound = nNodes * Omega * sqrt(density);
    
    % Interval size
    intervalSz = sqrt(2*epsilon / (eigenBound * nNodes));
    
    function qr = qRange(n)
        qr = A(n):intervalSz:(1 - B(n));
        if isempty(qr) % Whole interval is within intervalSz
            qr = [A(n) (1 - B(n))];
        elseif qr(end) ~= (1 - B(n)) % Did not include the endpoint
            qr = [qr (1 - B(n))]
        end
        assert(length(qr) >= 2, 'not enough intervals!');        
    end
            
    % Unaries: Term two of (Eq 4)
    D = cell(nNodes, 1);    
    for n = 1:nNodes
        qr = qRange(n);        
        D{n} = zeros(size(qr));                        
        
        t = theta(n);
        for iq = 1:length(qr)
            q = qr(iq);
            S = -q * log(q) - (1 - q) * log(1 - q);
            D{n}(iq) = -t * q + (deg(n) - 1) * S
        end        
    end
    
    % Pairwise terms (Eq 5)
    Vm = cell(nEdges, 1);
    vVec = zeros(nEdges, 1);
    for ne = 1:length(iVec)
        i = iVec(ne);
        j = jVec(ne);
        w = wVec(ne);
        vVec(ne) = ne;        
        
        alpha = exp(w) - 1;
        if alpha == 0
            warning('Got a zero weight; skipping. May cause problems.');            
            continue;
        end
                
        qir = qRange(i);
        qjr = qRange(j);
        Vm{ne} = zeros(length(qir), length(qjr));
        
        for iq = 1:length(qir)
            for jq = 1:length(qjr)
                q_i = qir(iq);
                q_j = qjr(jq);
                
                poly = [alpha -(1 + alpha*(q_i + q_j)) (1 + alpha)*q_i*q_j];
                rs = roots(poly);
                
                if alpha > 0
                    xi = min(rs);
                else 
                    xi = max(rs);
                end
                
                % (Eq 2); linearized
                marginal = [ 1 + xi - q_i - q_j
                             q_j - xi
                             q_i - xi
                             xi ];
                
                % Lemma 8: If neither q_i nor q_j are zero, then all
                % marginals should be positive.
                if q_i ~= 0 && q_j ~= 0
                    assert(all(marginal > 0));
                end
                S = -marginal' * log(marginal);
                
                Vm{ne}(iq, jq) = -w * xi - S;                
            end
        end        
    end
    
    % Symmetrize
    newWVec = ones(2 * nEdges, 1);
    
    siVec = vertcat(iVec, jVec);
    sjVec = vertcat(jVec, iVec);
    svVec = vertcat(vVec, vVec);
    newW = sparse(siVec, sjVec, newWVec, nNodes, nNodes);        
    Vi   = sparse(siVec, sjVec, svVec, nNodes, nNodes);        
end

