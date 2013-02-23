function V = makeMonge(N, M, figs)
% makeMonge Generate submodular matrix.
%
%   V = makeMonge(N, M, sigFigs) makes a submodular matrix of size NxM with
%   figs digits after the decimal place, if specified.
%
%   Currently it generates two random monotonic sequences and takes their
%   L1 norm. TODO: Vary functions.

    % Start by constructing a square matrix for simplicity
    maxL = max(N, M);    
    while true
        % Monotonic sequence
        %     Vm(:,:,2) = Vm(:,:,1).^2; % L2 term
        seq1 = cumsum(rand(1, maxL));
        seq2 = cumsum(rand(1, maxL));
        V = abs(bsxfun(@minus, seq1, seq2'));
        if nargin > 2 && figs
            V = fround(V, figs);
        end
        if IsMonge(V)            
            break;
        end        
    end
    
    % Slice to fit
    V = V(1:N,1:M);
    
    assert(IsMonge(V));
end
