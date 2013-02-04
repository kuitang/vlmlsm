function [ x, max_flow ] = MultiLabelSubModularBasic( D, W, V )
    
    [L, N] = size(D);
    nV = size(V, 3);
    E = nnz(W);
    
    % Edges for the graph
    % Nodes are zero-indexed
    nnodes = N*(L-1); % excluding s/t
    nedges = N*(L-2) + E*(L-1)*(L-1) + 2*nnodes;
    
    nodeS = nnodes;
    nodeT = nnodes + 1;
    
    ce = 1; % current edge
    ivec = zeros(nedges, 1);
    jvec = zeros(nedges, 1);
    ijvec = zeros(nedges, 1);
    jivec = zeros(nedges, 1);
            
    % Unary terms, source/sink weights, and stage ladder-climbing terms
    edgeCounter = 0;
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
            
            this_node = (r - 1) * (L - 1) + (k - 1);
            if qrk > 0     
                ivec(ce) = nodeS;
                jvec(ce) = this_node;               
                ijvec(ce) = qrk;                                
            else
                ivec(ce) = this_node;
                jvec(ce) = nodeT;                
                ijvec(ce) = -qrk;
            end
            assert(ijvec(ce) >= 0);
            ce = ce + 1;
                        
            if k <= L - 2 % between-states edge
                ivec(ce) = this_node;
                jvec(ce) = this_node + 1;
                ijvec(ce) = 0;
                jivec(ce) = 1e100; % Inf, but not literally
                
                ce = ce + 1;
                edgeCounter = edgeCounter + 1;
            end                                    
        end
    end
    assert(edgeCounter == N*(L - 2));
    edgeCounter = 0;
    
    % Pairwise terms
    for r = 1:N
        for rr = find(W(:,r))'
            fw = full(W(rr,r));
            w = real(fw);
            v = int32(imag(fw));
            assert(w > 0);
            assert(v > 0 && v <= nV);
            
            for k = 1:(L - 1)
                for kk = 1:(L - 1)
                    arr = w * (V(k,kk,v) ...
                               + V(k+1,kk+1,v) ...
                               - V(k+1,kk,v) ...
                               - V(k,kk+1,v));
                    arr = -arr / 2;
                    if abs(arr) <= eps
                        arr = 0;
                    end
                    assert(arr >= 0);
                    this_node = (r - 1) *  (L - 1) + (k - 1);
                    that_node = (rr - 1) * (L - 1) + (kk - 1);
                    ivec(ce) = this_node;
                    jvec(ce) = that_node;
                    ijvec(ce) = arr;
                    jivec(ce) = arr;
                    ce = ce + 1;                                        
                    edgeCounter = edgeCounter + 1;
                end
            end           
        end
    end
    assert(edgeCounter == E*(L-1)*(L-1));
    
    keep = ivec ~= jvec;
    ivec = ivec(keep);
    jvec = jvec(keep); 
    ijvec = ijvec(keep);
    jivec = jivec(keep);
    
    [max_flow, cut] = BK_mex(ivec, jvec, ijvec, jivec, nnodes);    
    
    % Translated
    x = zeros(1, N);
    for r = 1:N
        x(r) = L;
        breakout = false;
        for k = 1:L            
            if ~breakout
                % 1-index
                node = (r - 1) * (L - 1) + (k - 1) + 1;
                if ~cut(node) % node is in sink, not source
                    x(r) = k + 1;
                    breakout = true;
                end
            end
        end
    end
end

