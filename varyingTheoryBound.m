function [ dRange, theoryBound, a, b ] = varyingTheoryBound(rr)
% Gather theoretical bounds from varying R output
%
% [dRange, theoryBound] = varyingTheoryBound(rr)
%
% dRange      - Degree range
% theoryBound - Bounds calculated by getIntervalSz

    [A, B, ~]      = arrayfun(@(r) BBP(r.theta, r.W, 0, 0), rr, 'UniformOutput', false);
    [~, ~, bound, a, b] = cellfun(@(w, a, b) getIntervalSz(a, b, w, 0.01), {rr.W}, A, B);
            
    dMin = min([rr.d]);
    dMax = max([rr.d]);
    
    dRange = dMin:dMax;

    %theoryBound = arrayfun(@(e) bound(find([rr.d] == e, 1)), dMin:dMax);
    theoryBound  = arrayfun(@(e) mean(bound([rr.d] == e)), dRange);
    %arrayfun(@(e) var(bound([rr.d] == e)), dMin:dMax)
    %kappa       = arrayfun(@(e) kappa(find([rr.d] == e, 1)), dMin:dMax);
    %kappa        = arrayfun(@(e) mean(kappa([rr.d] == e)), dRange);
    %arrayfun(@(e) var(bound([rr.d] == e)), dMin:dMax)
    a            = arrayfun(@(e) mean(a([rr.d] == e)), dRange);
    b            = arrayfun(@(e) mean(b([rr.d] == e)), dRange);
end
