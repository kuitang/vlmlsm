function [ D, newW, Vi, Vm ] = boundMRFNew(theta, W, gams)
% boundMRF Construct an MRF for the Bethe bound approximation, for one-dim
%   [ D, newW, Vi, Vm, gams ] = boundMRF(theta, W, gams)
%   Output an MRF with N nodes and, each having as many labels as
%   discretization intervals determined by A(n) : (1 - B(n) : gamma.
%
%   theta - N-vector of unary weights
%   W     - N x N symmetric sparse matrix of pairwise edge weights (Eq 1)
%   gams  - cell array of vectors of discretized points
%
%   D, W, Vi, Vm - Arguments to pass to MultiLabelSubModular

    % From the symmetric graph, extract an edge list where iVec(n) < jVec(n)
    nNodes = length(theta);
    [iVec, jVec, wVec] = findUT(W);    
    nEdges = length(iVec);
    
    deg = zeros(nNodes, 1);
    for n = 1:nNodes
        deg(n) = sum(W(:,n) ~= 0);        
    end
    assert(all(deg > 0), 'All nodes must be connected');           
    
    assert(iscell(gams) && (length(gams) == nNodes), ...
        'boundMRFNew:gams', 'gams must be a N-cell');        
    
    % Unaries: Term two of (Eq 4)
    %
    % Reconciling the notation: (Eq 1) writes minus signs in front of the
    % summations, but MultiLabelSubModular expects positive-signed terms.
    % For the unary terms, we fold -theta into D.    
    D = cell(nNodes, 1);    
    for n = 1:nNodes        
        D{n} = zeros(size(gams{n}));                        
        
        % Node n in the constructed graph takes L_n states, where L_n is the
        % number of discrete levels at which we evaluate the energy.        
        t = theta(n);
        for iq = 1:length(gams{n})
            q = gams{n}(iq);
            
            m = [q 1 - q];
            
            ent = zeros(1,2);
                        
            totalNeg = sum(m(m < 0));
            assert(abs(totalNeg) < 1e-5, 'Too much negativity!');            
            
            ent(m <= 0) = 0;
            nz = m > 0;
            ent(nz) = -m(nz) .* log(m(nz));
                        
            D{n}(iq) = -t * q + (deg(n) - 1) * sum(ent);
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
        
        aij = exp(w) - 1;
        if aij == 0
            warning('Got a zero weight; skipping. May cause problems.');            
            continue;
        end
                
%         qir = gams{i};
%         qjr = gams{j};
%         Vm{ne} = zeros(length(qir), length(qjr));
        
%        Vm{ne}  = makePotential_mex(w, gams{i}, gams{j});        
        Vm{ne} = makePotential(w, gams{i}, gams{j});
%        assertElementsAlmostEqual(Vm{ne}, otherVm);
        
        assert(all(~vec(isnan(Vm{ne}))));
        assert(sum(vec(abs(imag(Vm{ne})))) == 0, 'Got some imagination!');
        
%         for iq = 1:length(qir)
%             for jq = 1:length(qjr)
%                 q_i = qir(iq);
%                 q_j = qjr(jq);
%                                 
%                 % Just inlined the goddamn thing
%                 % Compute with Welling and Teh's footnote 1 on pp. 55
%                 if aij == 0 % No interaction; independent marginals
%                     xi = q_i * q_j;
%                 else
%                     beta = 1 / aij;
%                     R = beta + q_i + q_j;
%                     xi = 0.5*(R - sign(beta)*sqrt(R^2 - 4*(1 + beta)*q_i*q_j));
%                 end
%                 
%                 marginal = [ 1 + xi - q_i - q_j, q_j - xi; q_i - xi, xi ];                                
%                 marginal = marginal(:);
%                 % Lemma 8: If neither q_i nor q_j are zero, then all
%                 % marginals should be positive.
%                 if q_i ~= 0 && q_j ~= 0
%                     assert(all(marginal > 0));
%                 end
%                 S = -marginal' * log(marginal);
%                 
%                 Vm{ne}(iq, jq) = -w * xi - S;                
%             end
%         end
                
%        assertElementsAlmostEqual(Vm{ne}, checkPot);
    end
    
    % Symmetrize
    newWVec = ones(2 * nEdges, 1);
    
    siVec = vertcat(iVec, jVec);
    sjVec = vertcat(jVec, iVec);
    svVec = vertcat(vVec, vVec);
    newW = sparse(siVec, sjVec, newWVec, nNodes, nNodes);        
    Vi   = sparse(siVec, sjVec, svVec, nNodes, nNodes);        
end

