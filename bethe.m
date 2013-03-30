function [ F ] = bethe( theta, W, oneMargs, twoMargs )
% bethe - Compute Bethe free energy given marginals
%   F = bethe(W, oneMargs, twoMargs)

    [iVec, jVec, wVec] = findUT(W);
    
    % Remember, entropy is NEGATIVE p log p!
    
    posProbs = squeeze(twoMargs(2,2,:));
    linMarg  = twoMargs(:);    
    twoF     = -sum(wVec .* posProbs) + sum(linMarg .* log(linMarg));
    
    degMin     = full(sum(W > 0, 2)) - 1;
    linOneMarg = [oneMargs, 1 - oneMargs];
    oneS       = -sum(linOneMarg .* log(linOneMarg), 2);
    oneF       = -sum(theta .* oneMargs) + sum(degMin .* oneS);
    
    F = oneF + twoF;    
end

