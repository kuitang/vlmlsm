function [xMin eMin] = MultiLabelSubModularBruteForce(D, W, Vi, Vm)
% Discrete optimization of the following multi-label functional
%
% x = argmin \sum_i D_i(x_i) + \sum_ij w_ij V_m (x_i, x_j)
%
% for N-dimensional x, where dimension m has K_m distinct labels and all
% interaction matrices V_m are submodular.
%
% Usage:
%   [x e] = MultiLabelSubModularBruteForce(D, G, Vi, Vm)
%
% Inputs:
%   D   - Unary costs. Length-N cell array of length-K_m vectors
%   W   - (Sparse) pair-wise weights (w_ij) of size NxN. Lower triangular
%         discarded.
%   Vi  - Interaction index matrix of size NxN. Entry ij selects the
%         interaction matrix for pair ij. Same shape and sparsity as W.
%   Vm  - Interaction matrices. Length-M cell array of size K_m x K_m
%
% Outputs:
%   xMin - optimal assignment
%   eMin - energy of optimal assignment. Total energy in e(1), unary energy
%          in e(2), interaction energy in e(3).
%
    N = length(D);    
    
    eMin  = [Inf Inf Inf];
    xMin = 0;    

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

%     assert( all( vij>= 1 ) && all( vij <= size(Vm,3) ),...
%         'MultiLabelSubModular:Vi',...
%         'Vi entries must be between 1 and %d', size(Vm,3));

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

    % Probe how many states each variable has
    limits = zeros(1, N);
    for n = 1:N
        limits(n) = length(D{n});
        assert(limits(n) >= 2, 'Each node must have at least two states!')
    end
    
    % Iterate in lexicographic order    
    xs = enumerate(limits);
    for nn = 1:size(xs, 1)
        x = xs(nn,:);
           
        e(3) = 0;
        % Pairwise term
        for r = 1:N  
            rrs = find(W(r,:));
            for rr = rrs(rrs > r)            
                w = full(W(r,rr));                
                                
                if w > 0
                    v = full(Vi(r,rr));
                    e(3) = e(3) + Vm{v}(x(r),x(rr));                    
                end
                
            end
        end
        
        e(2) = 0;
        for r = 1:N
            e(2) = e(2) + D{r}(x(r));
        end
                
        e(1) = sum(e(2:3));         
        if e(1) < eMin(1)            
            eMin = e;
            xMin = x;
        end 
       
    end    
end

