function [ M, xi ] = marginalize( alpha, q_i, q_j )
% Compute p(X_i = q_i, X_j = q_j).
%
% alpha - (i,j) entry of the alpha matrix (Below Eq 2)
% q_i   - value of the lower-numbered variable
% q_j   - value of the higher-numbered variable
    
    % Compute with Welling and Teh's footnote 1 on pp. 55    
    beta = 1/alpha;
    R = beta + q_i + q_j;
    xi = 0.5*(R - sign(beta)*sqrt(R^2 - 4*(1 + beta)*q_i*q_j));        

    p = [alpha -(1 + alpha*(q_i + q_j)) (1 + alpha)*q_i*q_j];
    rs = roots(p);

    if alpha > 0
        xiRoot = min(rs);
    else 
        xiRoot = max(rs);
    end
    
    assertElementsAlmostEqual(xi, xiRoot);
    
    M = [ 1 + xi - q_i - q_j, q_j - xi;
          q_i - xi,           xi ];

end

