%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);

% This parameter is definitely the "more sensitive" one.
%
% Swallow nNodes from the environment.
nJPts = 1000;
maxJ  = 5;
js    = linspace(0, maxJ, nJPts);

fn = ['unif_gridsearch_' num2str(nNodes) '_nodes_' datestr(now, 0);];

lineSearchParams = struct('nNodes', nNodes, 'etaMin', -1, 'etaMax', 1, 'nPts', 100, ...
                          'margThresh', 0.1, 'betheThresh', 0.1, 'nSeqRnd', 100, 'maxIntervals', 1e5, 'epsilon', 1e-2);

matlabpool open 12                

allProblems = cell(nJPts, 1);
parfor i = 1:nJPts
    allProblems{i} = unifLineSearch(js(i), lineSearchParams);
end

save fn;

matlabpool close

