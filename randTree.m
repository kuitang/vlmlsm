function [ T ] = randTree( nNodes )
% T = randTree(nNodes)
%
% Create a random uniform spanning tree T of nNodes nodes with Wilson's
% algorithm. Simplified for unweighted complete graph.
%
% T is a sparse nNodes x nNodes adjacency matrix.
%
% References
% [1] David Bruce Wilson, "Generating Random Spanning Trees More Quickly
%     than the Cover Time" STOC '96 doi>10.1145/237814.237880

    r = randsample(nNodes, 1);
    parents = zeros(nNodes, 1);
    inTree = false(nNodes, 1);
    inTree(r) = true;
    for i = 1:nNodes
        u = i;
        while ~inTree(u)
            pop = 1:nNodes;
            pop(u) = [];
            parents(u) = randsample(pop, 1);
            u = parents(u);
        end
        u = i;
        while ~inTree(u)
            inTree(u) = true;
            u = parents(u);
        end
    end
    
    % the root does not have a parent
    parentsIdxs = find(parents);
    parentsNz   = parents(parentsIdxs);
    wvec        = true(nNodes - 1, 1);
    
    % Symmetrize and sparsify
    T = sparse(vertcat(parentsIdxs, parentsNz), ...
               vertcat(parentsNz, parentsIdxs), ...
               vertcat(wvec, wvec));
end

