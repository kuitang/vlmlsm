T      = [0 1 2];
nNodes = [8 16 32];

TVec = vec(repmat(T, 4, 1));
nVec = repmat(nNodes, 1, 3);

% Do not run in parallel; we don't have the memory
for i = 1:9
   t = TVec(i);
   n = nVec(i);

   worstSparse(dW, t, n, 0.1);
end

