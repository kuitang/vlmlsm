function [x, e, elMat] = SingleLabelSubModular(D, W, Vi, Vm)
% Discrete optimization of the following single-label functional
%
% x = argmin \sum_i D_i(x_i) + \sum_ij w_ij V_k (x_i, x_j)
%
% for N-dimensional x, over 2 discrete labels
% Provided that the energy is multilabel submodular.
%
% Usage:
%   [x, e, elMat] = MultiLabelSubModular(D, W, Vi, Vm)
%
% Inputs:
%   D   - Data term of size Nx2
%   W   - Pair-wise weights (w_ij) of size NxN
%   Vi  - Index k for each pair ij, selecting the pair-wise interaction V_k. 
%         Vi is of size NxN
%         note that for each non-zero entry in W there must be a
%         corresponding non-zero entry in Vi. Entries must be proper
%         indices into Vm, i.e., integer numbers in the range
%         [1..size(Vm,3)]
%   Vm  - A concatanation of matrices V_k of size 2x2xK
%
% Outputs:
%   x   - optimal assignment
%   e   - energy of optimal assignment
%   elMat - the matrix passed to BK_mex for inspection
%
% This implements the construction on page 151 of
% [1] Kolmogorov, V. and Zabin, R.,
%     "What energy functions can be minimized via graph cuts?",
%     Pattern Analysis and Machine Intelligence, IEEE Transactions on,
%     Feb. 2004
%
% Based on code by Bagon Shai (see MultiLabelSubModular.m)

% check Monge
assert( IsMonge(Vm), 'MultiLabelSubModular:submodularityV',...
    'Matrices Vm are not submodular');

% make sure W is non negative
assert( all( nonzeros(W) > 0 ), 'MultiLabelSubModular:submodularityW',...
    'weights w_ij must be non-negative');

assert(size(D, 2) == 2, 'only binary labels supported!');

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
wvec = vertcat(swij + 1i*svij, swij + 1i*svij);
W = sparse(ivec, jvec, wvec, size(W,1), size(W,2));

nNodes = size(D, 1);

% Each node has an i->j edge, and each original edge shows up once. 
nEdges = nNodes + nnz(W) / 2;

giVec = zeros(nEdges, 1); gjVec = zeros(nEdges, 1);
gijVec = zeros(nEdges, 1); gjiVec = zeros(nEdges, 1);
ce = 1;

sNode = nNodes + 1;
tNode = nNodes + 2;

function addST(node, weight)
    if weight >= 0        
        giVec(ce) = sNode;
        gjVec(ce) = node;
        gijVec(ce) = weight;
    else
        giVec(ce) = node;
        gjVec(ce) = tNode;
        gijVec(ce) = -weight;
    end
    ce = ce + 1;    
end

function addEdge(i, j, ij, ji)
    giVec(ce) = i;
    gjVec(ce) = j;
    gijVec(ce) = ij;
    gjiVec(ce) = ji;
    ce = ce + 1;
end

% Add s/t edges for each node
for n = 1:nNodes
    w = D(n,2) - D(n,1);
    addST(n, w);    
end

% We have n < nn. Iterate first on nn by column and then select rows by n.
% I believe this is faster according to how sparse arrays are handled.
for nn = 1:nNodes
    ns = find(W(:,nn))';
    ns = ns(ns < nn);
    for n = ns
        fw = full(W(n,nn));
        w = real(fw);
        v = int32(imag(fw));

        % See Table 2
        wNST   = Vm(2,1,v) - Vm(1,1,v);
        wNNST  = Vm(2,2,v) - Vm(2,1,v);
        wNToNN = Vm(1,2,v) + Vm(2,1,v) - Vm(1,1,v) - Vm(2,2,v);

        addST(n, wNST);
        addST(nn, wNNST);

        assert(wNToNN >= 0);
        addEdge(n, nn, wNToNN, 0);        
    end
end

% Slow and stupid
gW = sparse(giVec, gjVec, gijVec);
[giVec, gjVec, gijVec] = find(gW);
gjiVec = zeros(length(giVec), 1);
elMat = horzcat(giVec, gjVec, gijVec, gjiVec);

[e, x] = BK_mex(giVec, gjVec, gijVec, gjiVec, nNodes, sNode, tNode);
% One indexing
x = x(2:end) + 1;

if nargout > 1 % compute energy as well
    N = size(D,1);

    e(3) = Vm(sub2ind(size(Vm), x(swi), x(swj), svij')) * swij; % pair-wise term
    e(2) = sum(D(sub2ind(size(D), 1:N, x))); % data term (unary)
    e(1) = sum(e(2:3));
end

end

%-------------------------------------------------------------------------%
function tf = IsMonge(M)
%
% Checks if matrix M is a Monge matrix.
%
% The following must hold for a Monge matrix (with finite entries):
%   M_ik + M_jl <= M_il + M_jk \forall 1 <= i <= j <= m, 1 <= k <= l <= n
%
% M may be 3D array, in which case each 2D "slice" is checked for Monge
% property
%

tf = true;

for si=1:size(M,3)
    [m n] = size(M(:,:,si));
    for ii=1:(m-1)
        for kk=1:(n-1)
            tf =  ( M(ii,kk) + M(ii+1,kk+1) <= M(ii, kk+1) + M(ii+1, kk) );
            if ~tf, return; end
        end
    end
end

end
