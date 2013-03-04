function [ theta, W ] = makeProblem( nNodes, a, b )
% makeProblem Make an MRF problem per Eq 1. W normalized to max 1 with
% unary weights drawn from a, b uniform.
%   [theta, W] = makeProblem(nNodes)

    theta = unifrnd(a, b, nNodes, 1);
    % W doesn't actually have to be Monge; it just needs to be positive
    W = rand(nNodes, nNodes);    
    W(1:nNodes+1:nNodes*nNodes) = 0;                     % Remove diagonal
    W = .5 * (W + W');                                   % Symmetrize    
    W = W ./ max(W(:));                                  % Normalize (VERY IMPORTANT TO PREVENT OVERFLOW!)
    W = sparse(W);

end

