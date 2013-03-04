function [x, e, elMat, maxFlow] = MultiLabelSubModular(D, W, Vi, Vm)
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
%   x - optimal assignment
%   e - energy of optimal assignment. Total energy in e(1), unary energy
%       in e(2), interaction energy in e(3).
%
%
% 
%-------------------------------------- 
%
% NOTE
% requires compilation of mex file
%
% >> mex -O -largeArrayDims -DNDEBUG graph.cpp maxflow.cpp...
%        BK_mex.cpp -output BK_mex
%
%--------------------------------------
%
% Copyright (c) Kui Tang
% Department of Applied Physics and Applied Mathematics
% Columbia University
%
% Based on code by Bagon Shai
%
% Copyright (c) Bagon Shai
% Department of Computer Science and Applied Mathmatics
% Wiezmann Institute of Science
% http://www.wisdom.weizmann.ac.il/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.
%
%
% If used in an academic research, the follwoing citation must be included
% in any resulting publication:
%
%   [1] D. Schlesinger and B. Flach, 
%       "Transforming an arbitrary MinSum problem into a binary one", 
%       Technical report TUD-FI06-01, Dresden University of Technology, April 2006.
%       http://www1.inf.tu-dresden.de/~ds24/publications/tr_kto2.pdf
%
%   [2] Yuri Boykov and Vladimir Kolmogorov,
%       "An Experimental Comparison of Min-Cut/Max-Flow Algorithms for 
%       Energy Minimization in Computer Vision",
%       PAMI, September 2004. 
% 
%   [3] Shai Bagon 
%       "Matlab Implementation of Schlezinger and Flach Submodular Optimization",
%       June 2012.
%

% check Monge
assert( IsMonge(Vm), 'MultiLabelSubModular:submodularityV',...
    'Matrices Vm are not submodular');

% make sure W is non negative
assert( all( nonzeros(W) > 0 ), 'MultiLabelSubModular:submodularityW',...
    'weights w_ij must be non-negative');

% remove diagonal of W
[swi, swj, swij] = findUT(W);
[svi, svj, svij] = findUT(Vi);

assert( isequal(svi, swi) && isequal(swj, svj), 'MultiLabelSubModular:Vi', ...
    'Index matrix Vi must have non-zeros pattern matching W');

assert( all(svij == round(svij)), 'MultiLabelSubModular:Vi',...
    'Vi must have integer index entries');

% assert( all( vij>= 1 ) && all( vij <= size(Vm,3) ),...
%     'MultiLabelSubModular:Vi',...
%     'Vi entries must be between 1 and %d', size(Vm,3));

% symmetrize: add the transpose entries
ivec = vertcat(swi, swj);
jvec = vertcat(swj, swi);
wvec = vertcat(swij + 1i*svij, swij + 1i*svij);
W = sparse(ivec, jvec, wvec, size(W,1), size(W,2));

% run the optimization
%xMex = MultiLabelSubModular_mex(D, W, Vm);
[x, maxFlow, elMat] = MultiLabelSubModularBasic(D, W, Vm);
%assert(all(xMex == x));

if nargout > 1 % compute energy as well
    N = length(D);    
        
    % Pairwise term
    e = zeros(1, 3);
    e(3) = 0;
    % Pairwise term
    for r = 1:N  
        rrs = find(W(r,:));
        for rr = rrs(rrs > r)            
            fw = full(W(r,rr));
            w = real(fw);
            v = int32(imag(fw));

            if w > 0
                v = full(Vi(r,rr));
                e(3) = e(3) + w * Vm{v}(x(r),x(rr));
            end

        end
    end
    
    e(2) = 0;
    for r = 1:N
        e(2) = e(2) + D{r}(x(r));
    end    
    
    e(1) = sum(e(2:3));
end

%-------------------------------------------------------------------------%
function tf = IsMonge(M)
%
% Checks if all entries of cell array M are Monge matrices.
%
% The following must hold for a Monge matrix (with finite entries):
%   M_ik + M_jl <= M_il + M_jk \forall 1 <= i <= j <= m, 1 <= k <= l <= n
%

tf = true;
TEN_EPS = 2.2204e-15;

N = length(M);
for n = 1:N
    V = M{n};
    [maxRow, maxCol] = size(V);
    maxRow = maxRow - 1;
    maxCol = maxCol - 1;
    for ii=1:maxRow
        for kk=1:maxCol
            diff = V(ii,kk) + V(ii+1,kk+1) - V(ii, kk+1) - V(ii+1, kk);
            if diff - 5 * TEN_EPS > 0
                tf = false;
                return
            end            
        end
    end
end
