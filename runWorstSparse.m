matlabpool open 12

T      = [0 2 4];
nNodes = [8 16 32];

TVec = vec(repmat(T, 4, 1));
nVec = repmat(nNodes, 1, 3);

parfor i = 1:12 
   t = TVec(i);
   n = nVec(i);

   fn = sprintf('bregular%d.%d.mat', n, log(n)/log(2));
   open(fn)

   worstSparse(dW, t, n, 0.1);
end

