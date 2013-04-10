function [ nearDisc ] = searchMRF( theta, W, A, B, intervalSz, trueMargs )
% boundMRF Construct an MRF for the Bethe bound approximation.
%   [ D, newW, Vi, Vm, qr ] = boundMRF( theta, W, A, B, alpha, epsilon )
%   Output an MRF with N nodes and, each having as many labels as
%   discretization intervals determined by A(n) : (1 - B(n) : gamma.
%
%   theta - N-vector of unary weights
%   W     - N x N symmetric sparse matrix of pairwise edge weights (Eq 1)
%   A     - N-vector of lower bounds
%   B     - N-vector of complementary upper bounds
%   alpha - output from BBP
%   epsilon - tolerance from the true optimum value
%
%   D, W, Vi, Vm - Arguments to pass to MultiLabelSubModular
%   qr           - N-cell array of probabilities corresponding the label levels
%
%   Equations and numbers taken from the 31 Dec 2012 draft.

    % From the symmetric graph, extract an edge list where iVec(n) < jVec(n)
    nNodes = length(theta);
    [iVec, jVec, wVec] = findUT(W);    
    nEdges = length(iVec);
    
    deg = zeros(nNodes, 1);
    for n = 1:nNodes
        deg(n) = sum(W(:,n) > 0);
        assert(deg(n) > 0, 'All nodes must be connected');
    end
            
    qr = cell(nNodes, 1);
    for n = 1:nNodes
        q = A(n):intervalSz:(1 - B(n));
        if isempty(q) % Whole interval is within intervalSz
            q = [A(n) (1 - B(n))];
        elseif q(end) ~= (1 - B(n)) % Did not include the endpoint
            q = [q (1 - B(n))];
        end
        
        assert(length(q) >= 2, 'not enough intervals!');                
        qr{n} = q;
    end
    
    % Search for closest marginals in the discretization
    nearDisc = zeros(size(trueMargs));
    for n = 1:nNodes
        margDist = abs(trueMargs(n) - qr{n});        
        [~, iMin] = min(margDist);
        nearDisc(n) = qr{n}(iMin);
    end        
       
end

