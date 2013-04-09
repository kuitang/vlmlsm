%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);

fn = ['rand_gridsearch_' num2str(nNodes) '_nodes_' datestr(now, 0);];

% Double-cycle star graph
% adj = [ 0 1 0 0 1 
%         1 0 1 1 1
%         0 1 0 1 0
%         0 1 0 1 0
%         0 1 1 0 1
%         1 1 0 1 0 ];

% Fully-connected nNodes graph
adj = ones(nNodes,nNodes) - diag(diag(ones(nNodes,nNodes)));

randSearchParams = struct('nNodes', nNodes, 'etaMin', -1, 'etaMax', 1, ...
                          'jMax', 2, 'plot', false, 'nIters', nIters, ... 
                          'adj', adj, 'margThresh', 0.1, 'epsilon', 0.1, ...
                          'betheThresh', 0.1, 'nSeqRnd', 100, 'maxIntervals', 1e5);

%% Or actually run it
matlabpool open 12
spmd
    %% One iter
    [allProblems, allFailVecs] = randSearch(randSearchParams);    
end

%save(fn);

matlabpool close
