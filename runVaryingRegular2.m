% RUN ON A SEPARATE MACHINE FOR EACH NNODES
% set dw so that the problems take nontrivial time

dw = 6;
outfn = sprintf('VaryingRegular_%d_RESULTS.mat', nNodes);
infn = sprintf('varyingregular%d.mat', nNodes);
load(infn);

matlabpool open 12
results = parVaryingRegular(graphs, dw);

save(outfn);

