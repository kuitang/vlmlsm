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
    
    % Iterate in lexicographic order
    xs = enumerate(L * ones(1, N));
    for nn = 1:size(xs, 1)
        x = xs(nn,:);
                
        % Pairwise term
        energy(3) = Vm(sub2ind(size(Vm), x(swi), x(swj), svij')) * swij;        
        energy(2) = sum(D(sub2ind(size(D), 1:N, x))); % data term (unary)        
        energy(1) = sum(energy(2:3));
        
        if energy(1) < emin(1)            
            emin = energy;
            xmin = x;
        end                        
    end    
end

