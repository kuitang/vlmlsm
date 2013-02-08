function [xmin emin] = MultiLabelSubModularBruteForce(D, W, Vi, Vm)
% Discrete optimization of the following multi-label functional
%
% x = argmin \sum_i D_i(x_i) + \sum_ij w_ij V_k (x_i, x_j)
%
% for N-dimensional x, over L discrete labels
% Provided that the energy is multilabel submodular.
%
% Usage:
%   [x e] = MultiLabelSubModularBruteForce(D, G, Vi, Vm)
%
% Inputs:
%   D   - Data term of size NxL
%   W   - Pair-wise weights (w_ij) of size NxN
%   Vi  - Index k for each pair ij, selecting the pair-wise interaction V_k. 
%         Vi is of size NxN
%         note that for each non-zero entry in W there must be a
%         corresponding non-zero entry in Vi. Entries must be proper
%         indices into Vm, i.e., integer numbers in the range
%         [1..size(Vm,3)]
%   Vm  - A concatanation of matrices V_k of size LxLxK
%
% Outputs:
%   x   - optimal assignment
%   e   - energy of optimal assignment
    [N, L] = size(D);
    
    emin  = [Inf Inf Inf];
    xmin = 0;    

    % remove diagonal of W
    [wi wj wij] = find(W);
    sel = wi < wj;
    if ~any(sel)
        sel = wi > wj;
    end
    [vi vj vij] = find(Vi);
    assert( isequal(vi, wi) && isequal(wj, vj), 'MultiLabelSubModular:Vi', ...
        'Index matrix Vi must have non-zeros pattern matching W');

    assert( all(vij == round(vij)), 'MultiLabelSubModular:Vi',...
        'Vi must have integer index entries');

    assert( all( vij>= 1 ) && all( vij <= size(Vm,3) ),...
        'MultiLabelSubModular:Vi',...
        'Vi entries must be between 1 and %d', size(Vm,3));

    swi = wi(sel);
    swj = wj(sel);
    swij = wij(sel);
    svij = vij(sel);

    % symmetrize: add the transpose entries
    ivec = vertcat(swi, swj);
    jvec = vertcat(swj, swi);
    wvec = vertcat(swij, swij);
    vvec = vertcat(svij, svij);
    W =  sparse(ivec, jvec, wvec, size(W,1), size(W,2));
    Vi = sparse(ivec, jvec, vvec, size(W,1), size(W,2));

    % Iterate in lexicographic order
    xs = enumerate(L * ones(1, N));
    for nn = 1:size(xs, 1)
        x = xs(nn,:);
           
        e(3) = 0;
        % Pairwise term
        for r = 1:N            
            for rr = (r+1):N                
                w = full(W(r,rr));                
                                
                if w > 0
                    v = full(Vi(r,rr));
                    e(3) = e(3) + Vm(x(r),x(rr),v);
                end
                
            end
        end
        
        e(2) = sum(D(sub2ind(size(D), 1:N, x))); % data term (unary)
        e(1) = sum(e(2:3));         
        if e(1) < emin(1)            
            emin = e;
            xmin = x;
        end 
       
    end    
end

