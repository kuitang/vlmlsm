% RUN ON A SEPARATE MACHINE FOR EACH NNODES
% set dw so that the problems take nontrivial time

w = 2;
outfn = sprintf('VaryingRegular_PTAS_%d_RESULTS.mat', nNodes);
infn = sprintf('varyingregular%d.mat', nNodes);
load(infn);

matlabpool open 12
%results = parVaryingRegularPTAS(graphs, dw);
% Use fixed w instead
results = parVaryingRegularPTAS3(graphs, w);


save(outfn);

