clear D W Vi Vm;


N = 1000;
abortOnError = true;

success = 0;
for n = 1:N

    [D, W, Vi, Vm] = gen_mrf(randi(4)+1, randi(4)+1, rand);

    % Ignore unaries
    %D = zeros(6, 3);

    % Symmetrize W
    %W = W - tril(W);
    %Vi = Vi - tril(Vi);

    [xb, eb] = MultiLabelSubModularBruteForce(D, W, Vi, Vm);
    [x, e, elMat] = MultiLabelSubModular(D, W, Vi, Vm);
            
    if abortOnError
        assert(all(x == xb));
        assertElementsAlmostEqual(eb, e);
    end
    
    if all(x == xb)
        success = success + 1;
    end
end

prob = success / N

% end
