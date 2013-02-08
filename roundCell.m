function [ Dc, Vmc ] = roundCell( p, D, Vm )
% p significant figures
    N = length(D);
    Dc = cell(N, 1);
    for n = 1:N
        Dc{n} = fround(D{n}, p);
    end

    M = length(Vm);
    Vmc = cell(M, 1);
    for m = 1:M
        Vmc{m} = fround(Vm{m}, p);
    end

end

