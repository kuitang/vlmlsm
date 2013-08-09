function [ M, xi ] = marginalize( alpha, q_i, q_j )
% [M, xi] = marginalize(alpha, q_i, q_j) Compute p(X_i = q_i, X_j = q_j).
%
%   alpha - (i,j) entry of the alpha matrix (Below Eq 2)
%   q_i   - value of the lower-numbered variable
%   q_j   - value of the higher-numbered variable
%
%   References
%   [1] Welling, Max and Teh, Yee Whye. "Belief optimization for binary
%       networks: a stable alternative to loopy belief propagation."
%       UAI 2001.
    
    fudge = 1e-10;

    % Compute with Welling and Teh's footnote 1 on pp. 55    
    if alpha == 0 % No interaction; independent marginals
        xi = q_i * q_j;
    else
        beta = 1/alpha;
        R = beta + q_i + q_j;
        
        discrim = R^2 - 4*(1 + beta)*q_i*q_j;
        if discrim < 0 && discrim > -fudge
            discrim = 0;
        end
        
        assert(discrim >= 0, 'Uh oh');
        
        xi = 0.5*(R - sign(beta)*sqrt(discrim));        
    end
        
    M = [ 1 + xi - q_i - q_j ; q_j - xi ; q_i - xi ; xi ];
end

