clear D W Vi Vm;

N = 1000;

for n = 1:N

    [D, W, Vi, Vm] = gen_mrf(randi(4)+1, randi(4)+1, rand);
    %[D, W, Vi, Vm] = gen_mrf(3, 2, rand);

    % Ignore unaries
    %D = zeros(6, 3);

    % Symmetrize W
    %W = W - tril(W);
    %Vi = Vi - tril(Vi);

%     [xb, eb] = MultiLabelSubModularBruteForce(D, W, Vi, Vm);
%     [x, e, elMat] = MultiLabelSubModular(D, W, Vi, Vm);
    %W = Vi > 0;
    [xb, eb] = MultiLabelSubModularBruteForce(D, W, Vi, Vm);
    [x, e, elMat] = MultiLabelSubModular(D, W, Vi, Vm);
    

        
    % There may be multiple unique solutions!
    assertElementsAlmostEqual(eb(1), e(1));
    if ~all(x == xb)
        disp('Non unique solution:');
        disp(['Brute force x = ' num2str(xb) ' energy = ' num2str(eb)]);
        disp(['Actual      x = ' num2str(x)  ' energy = ' num2str(e)]);
    end
    %assert(all(x == xb));
    %assertElementsAlmostEqual(eb, e);
end

