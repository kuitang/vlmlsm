%% Steup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);

%% Configure and set the problem
epsilon = .01;
theta = [1 2 3]';
nNodes = length(theta);
W = [0 1 2 ; 1 0 1 ; 2 1 0];
nEdges = nnz(W) / 2;

%% Calculate some brute forces
% Find the optimal of the Bethe free energy for the graph and ensure that
% our algorithm is within epsilon.

% theta, are x column vectors

% TODO: Generalize. (Something about squeeze; permute, {:} indexing)
E = @(x) -theta'*x - x'*triu(W)*x;
energies = zeros(2,2,2);
for xp = enumerate([2 2 2])'   
    x = xp - 1;
    energies(xp(1),xp(2),xp(3)) = E(x);
end
Ztrue = sum(exp(-energies(:)));
probs = exp(-energies) ./ Ztrue;

% Marginalize
bf12Marginal = squeeze(sum(probs, 3));
bf13Marginal = squeeze(sum(probs, 2));
bf23Marginal = squeeze(sum(probs, 1));

bfMarginal = zeros(nNodes, 1);
bfm = sum(bf13Marginal, 1);
bfMarginal(3) = bfm(2);
bfm = sum(bf12Marginal, 1);
bfMarginal(2) = bfm(2);
bfm = sum(bf23Marginal, 1);
bfMarginal(1) = bfm(2);

assertElementsAlmostEqual(sum(bf12Marginal, 2), sum(bf13Marginal, 2));
assertElementsAlmostEqual(sum(bf13Marginal, 1), sum(bf23Marginal, 1));
assertElementsAlmostEqual(sum(bf12Marginal, 1)', sum(bf23Marginal, 2));

%% Run our algorithm
[A, B, alpha] = BBP(theta, W);
[D, WW, Vi, Vm, qr] = boundMRF(theta, W, A, B, alpha, epsilon);
[x, e, elMat] = MultiLabelSubModular(D, WW, Vi, Vm);

% Translate levels back to marginals
oneMarginals = zeros(nNodes, 1);
for n = 1:nNodes
    oneMarginals(n) = qr{n}(x(n));
end

[siVec, sjVec, ~] = findUT(W);
nOutEdges = length(siVec);
assert(nOutEdges == nEdges);
twoMarginals = zeros(2, 2, nEdges);
for ne = 1:nEdges
    i = siVec(ne);
    j = sjVec(ne);
    
    a = alpha(i,j);
    q_i = oneMarginals(i);
    q_j = oneMarginals(j);
    
    twoMarginals(:,:,ne) = marginalize(a, q_i, q_j);
end

%% Test consistency of pairwise marginals
for ne = 1:nEdges
    i = siVec(ne);
    j = sjVec(ne);
    
    checkIMarginal = sum(twoMarginals(:,:,ne), 2);   
    checkJMarginal = sum(twoMarginals(:,:,ne), 1);
    
    assertElementsAlmostEqual(oneMarginals(i), checkIMarginal(2));
    assertElementsAlmostEqual(oneMarginals(j), checkJMarginal(2));
end
