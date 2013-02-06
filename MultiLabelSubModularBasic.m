function [ x, maxFlow, elMat ] = MultiLabelSubModularBasic( D, W, V )
    
    [L, N] = size(D);
    nV = size(V, 3);
    E = nnz(W); % Number of edges in the underlying graph.
        
    % Nodes are one-indexed
    nNodes = N*(L-1); % excluding s/t
    
    % This equation (except 2*nNodes) was lifted from the C++, but I'm not
    % sure if its correct, due to the alphas.
    %
    % Presumably, it's supposed to be the number of times you call add_edge,
    % excluding s and t edges.
    nEdges = N*(L-2) + E*(L-1)*(L-1) + nNodes;
    
    nodeDim = [N L-1];
    
    %sNode = nNodes + 1; tNode = nNodes + 2;
    sNode = -1 ; tNode = -2;

    ce = 1; % current edge
    iVec = zeros(nEdges, 1); jVec = zeros(nEdges, 1);
    ijVec = zeros(nEdges, 1); jiVec = zeros(nEdges, 1);
    
    function addST(node, weight)
        if weight >= 0        
            iVec(ce) = sNode;
            jVec(ce) = node;
            ijVec(ce) = weight;
        else
            iVec(ce) = node;
            jVec(ce) = tNode;
            ijVec(ce) = -weight;
        end
        ce = ce + 1;    
    end

    function addEdge(i, j, ij, ji)
        iVec(ce) = i;
        jVec(ce) = j;
        ijVec(ce) = ij;
        jiVec(ce) = ji;
        ce = ce + 1;
    end
    
    nZInfEdges = 0;    
    % Zero/infinity weights
    for r = 1:N
        for k = 1:(L - 2)
            lowNode  = sub2ind(nodeDim, r, k);
            highNode = sub2ind(nodeDim, r, k + 1);
            
            addEdge(lowNode, highNode, 0, 1e100);            
            nZInfEdges = nZInfEdges + 1;
        end
    end       
    assert(nZInfEdges == N * (L - 2));
    
    % Unary/source/sink weights
    nStEdges = 0;
    % For debugging
    allZeroQrkInner = true;    
    for r = 1:N
        for k = 1:(L - 1)
            qrk = 0;
            for rr = find(W(:,r))'
                fw = full(W(rr,r)); w = real(fw); v = int32(imag(fw));
                assert(w > 0);                
                assert(v > 0 && v <= nV);
                
                if r < rr
                    Vv = V(:,:,v);
                elseif r > rr
                    Vv = V(:,:,v)';
                else
                    assert('You cant have a self-edge');                    
                end
                
                qrk = qrk + w * (Vv(k,1) + Vv(k,L) - Vv(k+1,1) - Vv(k+1,L));
                %qrk = qrk + w * (Vv(k,L) - Vv(k+1,L));
            end

            qrk = qrk / 2;
            qrk = qrk + D(k,r) - D(k+1,r);
                        
            lowNode = sub2ind(nodeDim, r, k);
            addST(lowNode, qrk);            
            nStEdges = nStEdges + 1;            
        end
        if qrk ~= 0
            allZeroQrkInner = false;                
        end        
    end
    
    if allZeroQrkInner
        warning('All qrk inner sums were zeros. May have problems.');
    end
    
    assert(nStEdges == nNodes);    
    
    % Pairwise weights
    nAlphaEdges = 0;
    for r = 1:N
        for rr = find(W(:,r))'
            fw = full(W(rr,r)); w = real(fw); v = int32(imag(fw));
            assert(w > 0);
            assert(v > 0 && v <= nV);
            
            % By convention, k is the setting of r and kk is the setting of
            % rr. However, each V matrix-slice is indexed by convention
            % with lower-numbered variable on rows, higher-numbered
            % variable on columns. For instance, if V(:,:,v) contained the
            % x_1 and x_2 interactions and x_1 = k, x_2 = kk, then we would
            % look up V(k,kk,v).
            %
            % However, during the iteration, we will also encounter the
            % setting x_1 = kk, x_2 = k. Then we will need to look up
            % V(kk,k,v).
            %
            % Note that little-v is symmetric in r and rr.            
            if r < rr
                Vv = V(:,:,v);
            elseif r > rr
                Vv = V(:,:,v)';
            else
                assert('You cant have a self edge');
            end

            for k = 1:(L - 1)
                lowNode = sub2ind(nodeDim, r, k);                                                            
                
                for kk = 1:(L - 1)
                    highNode = sub2ind(nodeDim, rr, kk);
                                                                                    
                    arr = w * (Vv(k,kk) + Vv(k+1,kk+1) - Vv(k+1,kk) - Vv(k,kk+1));
                    %arr = -arr;
                    arr = -arr / 2;
                    % Sometimes we get minor negative perturbation
                    if arr < 0 && arr >= -eps
                        arr = 0;
                    end
                    assert(arr >= 0);
                    
                    addEdge(lowNode, highNode, arr, 0);    
                    %addEdge(lowNode, highNode, arr, arr);                    

                    nAlphaEdges = nAlphaEdges + 1;                    
                end
            end           
        end
    end    
    
    assert(all(iVec ~= jVec));    

    elMat = [iVec jVec ijVec jiVec];
    elMat(elMat == 1e100) = inf;
    
    [maxFlow, cut] = BK_mex(iVec, jVec, ijVec, jiVec, nNodes, sNode, tNode);
    cut = cut(2:end) % 1-index
    
    % Recall that x(n) = 1 if node n was assigned to the SINK.
    %
    % Let k' be the true optimal label. Then node (r, k) is labelled 1 (SINK)
    % if k >= k' and 0 otherwise. We assign k' as the first SINK label. See
    % the diagram on page 8.    
    x = zeros(1, N);
    for r = 1:N
        x(r) = L;        
        
        k = 1;
        % The highest value we explicitly track is L - 1.
        while k <= L - 1 && ~cut(sub2ind(nodeDim, r, k)) % while node is in SOURCE
            k = k + 1;
        end
        
        x(r) = k;        
    end
end

