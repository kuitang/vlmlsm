function [ A, B, alpha ] = BBP( theta, W, thresh, maxIter )
% BBP Bethe belief propagation (Algorithm 1)
%   [ A, B, alpha ] = BBP( theta, W, thresh, maxIter )
%
%   theta - Unary weights in original graph (Eq 1). COLUMN VECTOR
%   W     - (Sparse) symmetric edge weights in original graph
%
%   A     - Lower bounds
%   B     - Complementary upper bounds

    nNodes = size(W, 1);
    
    posW = zeros(nNodes, 1);
    negW = zeros(nNodes, 1);
    
    for j = 1:nNodes
        is = find(W(:,j));
        col = W(is,j);
        
        colSel = W(W(:,j) ~= 0,j);
        assert(all(col == colSel));
        posW(j) = sum(col(col > 0));
        negW(j) = -sum(col(col < 0));
    end
    
    A = sigmoid(theta - negW);
    B = 1 - sigmoid(theta + posW);        
    
    % I'M NOT SURE ABOUT THE alpha = exp(w) - 1 vs exp(abs(w))
    alpha = exp(abs(W)) - 1;
    
    converged = false;
    
    if nargin < 3
        thresh = .002;
    end
    if nargin < 4
        maxIter = Inf;
    end
    
    iter = 0;
    while ~converged && iter <= maxIter
        oldA = A;
        oldB = B;
        
        for i = 1:nNodes
            L = 1;
            U = 1;

            for j = find(W(i,:))
                a = alpha(i,j);
                if W(i,j) > 0
                    L = L * (1 + a*A(j) / (1 + a*(1 - B(i))*(1 - A(j))));
                    U = U * (1 + a*B(j) / (1 + a*(1 - A(i))*(1 - B(j))));
                else
                    L = L * (1 + a*B(j) / (1 + a*(1 - B(i))*(1 - B(j))));
                    U = U * (1 + a*A(j) / (1 + a*(1 - A(i))*(1 - A(j))));
                end
            end
                                     
            A(i) = 1 / (1 + exp(-theta(i) + negW(i)) / L);
            B(i) = 1 / (1 + exp(theta(i)  + posW(i)) / U);
            
            % Lemma 9: At every iteration, each element of A, B
            % monotonically increase   
            assert(A(i) > oldA(i));
            assert(B(i) > oldB(i));
        end
        
        dA = abs(A - oldA);
        dB = abs(B - oldB);
        if all(dA < thresh) && all(dB < thresh)
            converged = true;
        elseif iter >= maxIter
            warning(['Convergence threshold ' num2str(thresh) ...
                     ' not reached after ' num2str(maxIter) ' iterations']);            
        end
        
        iter = iter + 1;
    end
end

