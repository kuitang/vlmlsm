function [ x, maxFlow, elMat ] = MultiLabelSubModularBasic( D, W, V )
% x = MultiLabelSubModularBasic(D, W, V) core computation for
%     the MultiLabelSubModular routine.
%
% Inputs:
%   D   - Unary costs. Length-N cell array of length-K_m vectors
%   W   - (Sparse) complex symmetric weight and index matrix of size NxN.
%         Real part is the edge weight. Imaginary part is an integral index
%         into V, selecting the interaction matrix to use.
%   V   - Interaction matrices. Length-M cell array of size K_m x K_m
%
% Outputs:
%   x       - Solution vector
%   maxFlow - Maximum flow in the constructed graph
%   elMat   - Edge list for the constructed graph
    
    N = length(D);
    M = length(V);                
    
    E = nnz(W); % Number of edges in the underlying graph.
    
    % Convert (r,k) indices to linear indices. offset(r) is the linear
    % index of (r,1) MINUS ONE. To compute index of (r,k), write
    %
    % offset(r) + k.
    
    nStates = zeros(N, 1);
    for n = 1:N
        nStates(n) = length(D{n});
        assert(nStates(n) >= 2, 'Each node must have at least two states!')
    end
    maxNStates = max(nStates);
    nStatesMinusOne = nStates - 1;
    offset = cumsum(nStatesMinusOne); % shift it
    offset = [ 0; offset(1:(end-1)) ];

    
    % Each node r nStatesMinusOne(r) nodes in the graph.
    nNodes = sum(nStatesMinusOne); % excluding s/t    

    % Upper bound the number of edges.
    % TODO: Do something smarter if this runs into memory constraints. Or
    % maybe just rely on BK's dynamic allocation (trade time for memory)
    nEdges = nNodes + sum(nStates - 2) + E * max(nStates);
    
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
        assert(i ~= j, 'Cannot add self-loop');
        assert(i <= nNodes && j <= nNodes, 'Node index out of bounds');
        iVec(ce) = i;
        jVec(ce) = j;
        ijVec(ce) = ij;
        jiVec(ce) = ji;
        ce = ce + 1;
    end
    
    nZInfEdges = 0;    
    % Zero/infinity weights
    for r = 1:N
        for k = 1:(nStates(r) - 2)
            lowNode  = offset(r) + k;
            highNode = offset(r) + k + 1;            
            
            addEdge(lowNode, highNode, 0, 1e100);            
            nZInfEdges = nZInfEdges + 1;
        end
    end       
    assert(nZInfEdges == sum(nStates - 2));
    
    % Unary/source/sink weights
    nStEdges = 0;        
    for r = 1:N
        for k = 1:nStatesMinusOne(r)        
            qrk = 0;
            for rr = find(W(:,r))'
                fw = full(W(rr,r)); w = real(fw); v = int32(imag(fw));
                assert(w > 0);                
                assert(v > 0 && v <= M);
                
                % By convention, k is the setting of r and kk is the setting of
                % rr. Each V matrix is indexed by lower-numbered variable
                % on rows, higher-numbered variable on columns.
                %
                % However, during the iteration, we will also encounter the
                % setting x_1 = kk, x_2 = k. Then we will need to look up
                % V{v}(kk,k) instead.
                %
                % Note that little-v is symmetric in r and rr by definition
                % of the W matrix. But the matrix V{v} is note symmetric.
                
                if r < rr
                    Vv = V{v};
                elseif r > rr
                    Vv = V{v}';
                else
                    assert('You cant have a self-edge');                    
                end
                
                % Dimensions of the interaction matrix must match
                assert(nStates(r)  == size(Vv, 1));
                assert(nStates(rr) == size(Vv, 2));
                
                qrk = qrk + w * (Vv(k,1) + Vv(k,end) - Vv(k+1,1) - Vv(k+1,end));                
            end

            qrk = qrk / 2;
            qrk = qrk + D{r}(k) - D{r}(k+1);
            
            node = offset(r) + k;            
            addST(node, qrk);            
            nStEdges = nStEdges + 1;            
        end       
    end    
    
    assert(nStEdges == nNodes);    
    
    % Pairwise weights
    nAlphaEdges = 0;
    for r = 1:N
        for rr = find(W(:,r))'
            fw = full(W(rr,r)); w = real(fw); v = int32(imag(fw));
            assert(w > 0);
            assert(v > 0 && v <= M);
                        
            if r < rr
                Vv = V{v};
            elseif r > rr
                Vv = V{v}';
            else
                assert('You cant have a self edge');
            end

            for k = 1:nStatesMinusOne(r)
                lowNode = offset(r) + k;                
                
                for kk = 1:nStatesMinusOne(rr)
                    highNode = offset(rr) + kk;                    
                    
                    inside = Vv(k,kk) + Vv(k+1,kk+1) - Vv(k+1,kk) - Vv(k,kk+1);
                    if abs(inside) <= 10*eps
                        inside = 0;
                    end
                    arr = -(w * inside) / 2;                                        
                    assert(arr >= 0);
                    
                    addEdge(lowNode, highNode, arr, 0);                        
                    nAlphaEdges = nAlphaEdges + 1;                    
                end
            end           
        end
    end    
    
    assert(all(iVec(iVec ~= 0) ~= jVec(jVec ~= 0)));    

    elMat = [iVec jVec ijVec jiVec];
    elMat(elMat == 1e100) = inf;
    
    [maxFlow, cut] = BK_mex(iVec, jVec, ijVec, jiVec, nNodes, sNode, tNode);
    cut = cut(2:end); % 1-index
    
    % Recall that x(n) = 1 if node n was assigned to the SINK.
    %
    % Let k' be the true optimal label. Then node (r, k) is labelled 1 (SINK)
    % if k >= k' and 0 otherwise. We assign k' as the first SINK label. See
    % the diagram on page 8.    
    x = zeros(1, N);
    for r = 1:N
        x(r) = nStates(r);        
        
        k = 1;
        % The highest value we explicitly track is nStatesMinusOne(r)
        while k <= nStatesMinusOne(r) && ~cut(offset(r) + k) % while node is in SOURCE
            k = k + 1;
        end
        
        x(r) = k;        
    end
end
