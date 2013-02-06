clear D W Vi Vm;

success = 0;
% while true
    
    [D, W, Vi, Vm] = gen_mrf(3, 3, .5);

    % Ignore unaries
    %D = zeros(6, 3);

    % Symmetrize W
    %W = W - tril(W);
    %Vi = Vi - tril(Vi);

    [xb, eb] = MultiLabelSubModularBruteForce(D, W, Vi, Vm)
    [x, e, elMat]   = MultiLabelSubModular(D, W, Vi, Vm);
    
    % Eventually this will fail.
    assert(all(x == xb));
    assertElementsAlmostEqual(eb, e);
    
    success = success + 1

% end
