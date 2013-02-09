function [ M, xi ] = marginalize( alpha, q_i, q_j )
% Compute p(X_i = q_i, X_j = q_j).
%
% alpha - (i,j) entry of the alpha matrix (Below Eq 2)
% q_i   - value of the lower-numbered variable
% q_j   - value of the higher-numbered variable

    p = [alpha -(1 + alpha*(q_i + q_j)) (1 + alpha)*q_i*q_j];
    rs = roots(p);

    if alpha > 0
        xi = min(rs);
    else 
        xi = max(rs);
    end
    
    M = [ 1 + xi - q_i - q_j, q_j - xi;
          q_i - xi,           xi ];

end

