function [ gams ] = fillGams(gamma, A, B, oddMethod)
% fillGams  Return cell-array of schedule for a fixed (per-dimension) gamma

    if nargin < 4
        oddMethod = false;
    end

    N = length(A);
    if length(gamma) == 1
        gamma = gamma * ones(N, 1);
    end
    
    gams = cell(N, 1);
    
    for i = 1:N
        mb = 1 - B(i);                
        
        if oddMethod
            g = (A(i) + gamma(i)) : (2*gamma(i)) : mb;
        else
            g = A(i):gamma(i):mb;
        end        
        
        if length(g) == 0
            g = [A(i) mb];
        elseif g(end) ~= mb
            g(end+1) = mb;
        elseif length(g) == 1
            g(end+1) = g(end) + min(eps, gamma(i));
        end
        gams{i} = g;
    end

end
