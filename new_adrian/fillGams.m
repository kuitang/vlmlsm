function [ gams ] = fillGams(gamma, A, B)
% fillGams  Return cell-array of schedule for a fixed (per-dimension) gamma

    N = length(A);
    if length(gamma) == 1
        gamma = gamma * ones(N, 1);
    end
    
    gams = cell(N, 1);
    
    for i = 1:N
        mb = 1 - B(i);
        g = A(i):gamma(i):mb;
        if g(end) ~= mb
            g(end+1) = mb;
        end
        gams{i} = g;
    end

end
