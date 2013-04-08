%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);

% This parameter is definitely the "more sensitive" one.
%
% Swallow nNodes from the environment.
nJPts = 10000;
maxJ  = 10;
js    = linspace(eps, maxJ, nJPts);

fn = ['unif_gridsearch_' num2str(nNodes) '_nodes_' datestr(now, 0);];

% Double-cycle star graph
% adj = [ 0 1 0 0 1 
%         1 0 1 1 1
%         0 1 0 1 0
%         0 1 0 1 0
%         0 1 1 0 1
%         1 1 0 1 0 ];

adj = ones(nNodes,nNodes) - diag(diag(ones(nNodes,nNodes)));
lineSearchParams = struct('nNodes', size(adj,1), 'etaMin', -5, 'etaMax', 5, 'nPts', 10000, 'adj', adj, 'plot', false, ...
                          'margThresh', 0.1, 'betheThresh', 0.1, 'nSeqRnd', 100, 'maxIntervals', 1e5, 'epsilon', 1e-2);

%% Or actually run it
matlabpool open 12
allProblems = cell(nJPts, 1);
failVecs    = cell(nJPts, 1);
parfor i = 1:nJPts
    %% Just one!
    [allProblems{i}, failVecs{i}] = unifLineSearch(js(i), lineSearchParams);
end

save(fn);

matlabpool close

