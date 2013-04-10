function [ A ] = randGraph(nNodes, deg)
% randGraph - Generate a Renyi graph from one with average deg degrees
%
% A      - adjacency matrix
% nNodes - nodes
% deg    - average degree

    while true
        pEdge = deg / (nNodes - 1);
        A = rand(nNodes, nNodes) < pEdge;

        % Remove the diagonal
        A = A - diag(diag(A));
        % Keep only the upper triangle
        A = triu(A);
        A = sparse(A + A');

        [nComponents, ~] = graphconncomp(A);
        
        if nComponents == 1
            break
        end
    end    
    
end

