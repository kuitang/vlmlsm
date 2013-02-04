function [ D, W, Vi, Vm ] = gen_mrf( N, L, edge_prob )
% [D, W, Vi, Vm] = gen_mrf(N, L, Nmonge, edge_prob)
% Generate test-case submodular mrf for MultiLabelSubModular
% [D, W, Vi, Vm] inputs to MultiLablSubModular_mex
%
% N - number of nodes
% L - number of labels
% Nmonge - number of Monge matrices to create
% edge_prob of all N*(N-1)/2 edges, what fraction will appear?

    D = rand(N, L);
    
    nzmat = rand(N, N) < edge_prob;
    W = ones(N, N) .* nzmat;
    %W = abs(rand(N, N) .* nzmat);
  
    seq = cumsum(rand(1, L));
    Vm = abs(bsxfun(@minus, seq, seq'));
    Vi = ones(N, N);
    
%     Vm = zeros(L, L, 2);
%     % Generate a random monotonic sequence
%     seq = cumsum(rand(1, L));    
%     Vm(:,:,1) = abs( bsxfun(@minus, seq, seq') ); % L1 term
%     Vm(:,:,2) = Vm(:,:,1).^2; % L2 term
    assert(IsMonge(Vm));
    
    
    %Vi = randi(2, N, N) .* nzmat;
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
