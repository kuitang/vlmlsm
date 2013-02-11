%% Steup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 200;

totDiffs = zeros(nTrials, 1);
epsilon = 1e-5;

%%
for t = 1:nTrials
    theta = -2*rand(3, 1) + 1;    
    nNodes = length(theta);
    W = rand(3, 3);
    W(1:3+1:3*3) = 0;   
    W = .5 * (W + W');
    nEdges = nnz(W) / 2;
    assert(nEdges == 3);

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
    bfTwoMarginals = zeros(2,2,3);
    bfTwoMarginals(:,:,1) = squeeze(sum(probs, 3)); % 1,2
    bfTwoMarginals(:,:,2) = squeeze(sum(probs, 2)); % 1,3
    bfTwoMarginals(:,:,3) = squeeze(sum(probs, 1)); % 2,3

    bfMarginals = zeros(nNodes, 1);
    bfm = sum(bfTwoMarginals(:,:,1), 2); % 1
    bfMarginals(1) = bfm(2);
    bfm = sum(bfTwoMarginals(:,:,1), 1); % 2
    bfMarginals(2) = bfm(2);
    bfm = sum(bfTwoMarginals(:,:,2), 1); % 3
    bfMarginals(3) = bfm(2);

    % assertElementsAlmostEqual(sum(bf12Marginal, 2), sum(bf13Marginal, 2));
    % assertElementsAlmostEqual(sum(bf13Marginal, 1), sum(bf23Marginal, 1));
    % assertElementsAlmostEqual(sum(bf12Marginal, 1)', sum(bf23Marginal, 2));

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

    totDiffs(t) = sum(abs(bfMarginals - oneMarginals));
    
%     disp(['Optimization finished with epsilon = ' num2str(epsilon)]);
%     disp(['One-norm of one-marginals: ' num2str(sum(abs(bfMarginals - oneMarginals))) ]);
%     disp(['Two-norm of one-marginals: ' num2str(norm(bfMarginals - oneMarginals)) ]);
%     disp(['One-norm of two-marginals: ' num2str(sum(abs(bfTwoMarginals(:) - twoMarginals(:)))) ]);
%     disp(['Two-norm of two-marginals (lineralized): ' num2str(norm(bfTwoMarginals(:) - twoMarginals(:))) ]);

end

figure;
hist(totDiffs);
title(['1-norm differences of brute force and computed marginals for \epsilon = ' num2str(epsilon)]);
