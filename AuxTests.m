function [margs, margEnts, binEnts] = AuxTests(vals, alphas)
% Compare with AuxTests_mex.

    N = length(vals);
    assert(length(alphas) == N^2);
    margEnts = zeros(N^2, 1);    
    margs    = zeros(4, N^2);

    binEnts = -vals .* log(vals) - (1 - vals) .* log(1 - vals);
    
    ind = 1;
    for i = 1:N
        for j = 1:N            
            [marg, ~] = marginalize(alphas(ind), vals(i), vals(j));
            marg = marg(:);
            margs(:,ind) = marg;
            margEnts(ind) = -marg' * log(marg);
            
            ind = ind + 1;
        end
    end
    
    binEnts = binEnts';
end

