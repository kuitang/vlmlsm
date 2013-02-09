function [ D, newW, Vi, Vm ] = BoundMRF( theta, W, A, B, gamma )
% Construct an MRF for the Bethe bound approximation.
%
% Output an MRF with N nodes and, each having as many labels as
% discretization intervals determined by A(n) : (1 - B(n) : gamma.
%
% theta - N-vector of unary weights
% W     - N x N symmetric sparse matrix of pairwise edge weights (Eq 1)
% A     - N-vector of lower bounds
% B     - N-vector of complementary upper bounds
% gamma - interval size
%
% D, W, Vi, Vm - Arguments to pass to MultiLabelSubModular
%
% Equations and numbers taken from the 31 Dec 2012 draft.

    % From the symmetric graph, extract an edge list where iVec(n) < jVec(n)
    nNodes = find(W, 1);
    [iVec, jVec, wVec] = find(W);
    sel = iVec < jVec;
    iVec = iVec(sel);
    jVec = jVec(sel);
    wVec = wVec(sel);
    
    nEdges = length(iVec);
            
    % Unaries: Term two of (Eq 4)
    D = cell(nNodes, 1);    
    for n = 1:nNodes
        qRange = A(n):gamma:(1 - B(n));
        D{n} = zeros(size(qRange));
        
        deg = sum(W(n,:) > 0);
        assert(deg >= 1, 'Each node must have degree >= 1');
        
        for iq = 1:length(qRange)
            q = qRange(iq);
            S = -q * log(q) - (1 - q) * log(1 - q);
            D{n}(iq) = -theta * q + (deg - 1) * S
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
                
        qiRange = A(i):gamma:(1 - B(i));
        qjRange = A(j):gamma:(1 - B(j));
        Vm{ne} = zeros(length(qiRange), length(qjRange));
        
        for iq = 1:length(qiRange)
            for jq = 1:length(qjRange)
                q_i = qiRange(iq);
                q_j = qjRange(jq);
                
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
    newW = sparse(siVec, sjVec, newWVec, nEdges, nEdges);        
    Vi   = sparse(siVec, sjVec, svVec, nEdges, nEdges);        
end

