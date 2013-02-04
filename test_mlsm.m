[D, W, Vi, Vm] = gen_mrf(5, 2, 1);

% Ignore unaries
%D = zeros(6, 3);

% Symmetrize W
W = W - tril(W);
Vi = Vi - tril(Vi);

[xb, eb] = MultiLabelSubModularBruteForce(D, W, Vi, Vm)
[x, e]   = MultiLabelSubModular(D, W, Vi, Vm)
assert(all(x == xb));
assertElementsAlmostEqual(eb, e);
