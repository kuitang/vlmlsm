matlabpool open 12

T      = [0 2 4];
nNodes = [8 16 32 64];

TVec = vec(repmat(T, 4, 1));
nVec = repmat(nNodes, 1, 3);


allGraphs = {};
for n = nNodes
    fn = sprintf('bregular%d.%d.mat', n, log(n)/log(2));
    load(fn);
    allGraphs{n} = graphs;
end

parfor i = 1:12
   t = TVec(i);
   n = nVec(i);

   for dw = [1 2 4 8 16 32]
       randSparse(dw, t, allGraphs{n});
   end
end

