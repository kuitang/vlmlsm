%-------------------------------------------------------------------------%
function tf = IsMonge(M)
%
% Checks if matrix M is a Monge matrix.
%
% The following must hold for a Monge matrix (with finite entries):
%   M_ik + M_jl <= M_il + M_jk \forall 1 <= i <= j <= m, 1 <= k <= l <= n
%
% M may NOT be a 3-D array.

tf = true;

[m n] = size(M);
for ii=1:(m-1)
    for kk=1:(n-1)
        tf =  ( M(ii,kk) + M(ii+1,kk+1) <= M(ii, kk+1) + M(ii+1, kk) );
        if ~tf, return; end
    end
end
