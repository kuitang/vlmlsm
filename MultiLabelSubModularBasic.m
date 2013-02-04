function [ x, max_flow ] = MultiLabelSubModularBasic( D, W, V )
    
    [L, N] = size(D);
    nV = size(V, 3);
    E = nnz(W);
    
    % Edges for the graph
    % Nodes are one-indexed
    nNodes = N*(L-1); % excluding s/t
    nEdges = N*(L-2) + E*(L-1)*(L-1) + 2*nNodes;
    
    nodeDim = [N L-1];
    
    sNode = nNodes + 1;
    tNode = nNodes + 2;
    
    ce = 1; % current edge
    iVec = zeros(nEdges, 1);
    jVec = zeros(nEdges, 1);
    ijVec = zeros(nEdges, 1);
    jiVec = zeros(nEdges, 1);
            
    % Unary terms, source/sink weights, and stage ladder-climbing terms    
    nStEdges = 0;
    nZInfEdges = 0;
    nAlphaEdges = 0;
    for r = 1:N
        for k = 1:(L - 1)
            qrk = 0;
            for rr = find(W(:,r))'
                fw = full(W(rr,r));
                w = real(fw);
                v = int32(imag(fw));
                assert(w > 0);                
                assert(v > 0 && v <= nV);

                qrk = qrk + w * (V(k,1,v) ...
                                 + V(k,L,v) ...
                                 - V(k+1,1,v) ...
                                 - V(k+1,L,v));  
            end
            qrk = qrk / 2;
            qrk = qrk + D(k,r) - D(k+1,r);
                        
            this_node = sub2ind(nodeDim, r, k);            
            if qrk > 0     
                iVec(ce) = sNode;
                jVec(ce) = this_node;               
                ijVec(ce) = qrk;                                
            else
                iVec(ce) = this_node;
                jVec(ce) = tNode;                
                ijVec(ce) = -qrk;
            end
            assert(ijVec(ce) >= 0);
            nStEdges = nStEdges + 1;
            ce = ce + 1;
                        
            if k <= L - 2 % between-states edge
                iVec(ce) = this_node;
                jVec(ce) = this_node + 1;
                ijVec(ce) = 0;
                jiVec(ce) = 1e100; % Inf, but not literally
                
                ce = ce + 1;
                nZInfEdges = nZInfEdges + 1;                
            end                                    
        end
    end    
    
    % Pairwise terms
    for r = 1:N
        for rr = find(W(:,r))'
            fw = full(W(rr,r));
            w = real(fw);
            v = int32(imag(fw));
            assert(w > 0);
            assert(v > 0 && v <= nV);
            
            for k = 1:(L - 1)
                this_node = sub2ind(nodeDim, r, k);                                                            
                for kk = 1:(L - 1)
                    arr = w * (V(k,kk,v) ...
                               + V(k+1,kk+1,v) ...
                               - V(k+1,kk,v) ...
                               - V(k,kk+1,v));
                    arr = -arr / 2;
                    if arr < 0 && arr >= -eps
                        arr = 0;
                    end
                    assert(arr >= 0);                    
                    
                    that_node = sub2ind(nodeDim, rr, kk);                    
                    iVec(ce) = this_node;
                    jVec(ce) = that_node;
                    ijVec(ce) = arr;
                    jiVec(ce) = arr;
                    ce = ce + 1;     
                    nAlphaEdges = nAlphaEdges + 1;                    
                end
            end           
        end
    end    
    
    keep = iVec ~= jVec;
    iVec = iVec(keep);
    jVec = jVec(keep); 
    ijVec = ijVec(keep);
    jiVec = jiVec(keep);
    
    % Translate to zero-index
    [max_flow, cut] = BK_mex(iVec - 1, jVec - 1, ijVec, jiVec, nNodes);    
    
    cut
    % Translated
    x = zeros(1, N);
    for r = 1:N
        x(r) = L;
        breakout = false;
        for k = 1:(L-1)
            if ~breakout
                % 1-index
                node = sub2ind(nodeDim, r, k);                
                if cut(node) % node is in sink, not source
                    x(r) = k;
                    breakout = true;
                end
            end
        end
    end
end

