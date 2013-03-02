function [ theta, W ] = makeProblem( nNodes )
% makeProblem Make an MRF problem per Eq 1. W normalized to max 1.
%   [theta, W] = makeProblem(nNodes)

    theta = -nNodes/2*rand(nNodes, 1);
    W = makeMonge(nNodes, nNodes);    
    W(1:nNodes+1:nNodes*nNodes) = 0;                     % Remove diagonal
    W = .5 * (W + W');                                   % Symmetrize    
    W = W ./ max(W(:));                                  % Normalize (VERY IMPORTANT TO PREVENT OVERFLOW!)
    W = sparse(W);

end

