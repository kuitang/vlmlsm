theta = [1 2 3]';
W = [0 1 2 ; 1 0 1 ; 2 1 0];

% theta, are x column vectors
E = @(x) -theta'*x - x'*triu(W)*x;

Ztrue = 0;
for xp = enumerate([2 2 2])'
    x = xp - 1;
    Ztrue = Ztrue + exp(-E(x));
end

[A, B, alpha] = BBP(theta, W);
epsilon = .1;

[D, WW, Vi, Vm] = boundMRF(theta, W, A, B, alpha, epsilon);

[x, e, elMat] = MultiLabelSubModular(D, WW, Vi, Vm);
