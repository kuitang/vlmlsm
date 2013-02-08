function [ D, W, Vi, Vm ] = gen_mrf( N, maxL, edge_prob, sigFigs )
% [D, W, Vi, Vm] = gen_mrf(N, L, Nmonge, edge_prob)
% Generate test-case submodular mrf for MultiLabelSubModular
%
% N         - number of nodes
% maxL      - Maximum number of labels
% edge_prob - of all N*(N-1)/2 edges, what fraction will appear?
% sigFigs   - if present, round (generate nice numbers for debug cases)

    if nargin <= 3
        sigFigs = 0; % Sentinel
    end
    
    D = cell(N, 1);
    for n = 1:N
        % Must have at least two states!
        L = 1 + randi(maxL - 1, 1);
        D{n} = rand(1, L);
        if sigFigs
            D{n} = fround(D{n}, sigFigs);
        end
    end    
    
    nzmat = 0;
    while sum(nzmat(:)) == 0
        nzmat = rand(N, N) < edge_prob;
        % remove diagonal;       
        % http://www.mathworks.com/matlabcentral/newsreader/view_thread/47182
        nzmat(1:N+1:N*N) = 0;   
    end
    
    W = ones(N, N) .* nzmat;        
    %W = abs(rand(N, N) .* nzmat);
    % Symmetrize
    W = 0.5 * (W + W');
  
    % For each nonzero in W (existing edge), create a Monge matrix of the
    % appropriate dimensionality    
        
    Vi = zeros(N, N);
    
    [i, j] = find(W);    
    nNZ = length(i);
    Vm = cell(nNZ, 1);
    
    iMonge = 1;
    for n = 1:nNZ
        ii = i(n);
        jj = j(n);
        if ii < jj
            Li = length(D{ii});
            Lj = length(D{jj});
            
            Vm{iMonge} = makeMonge(Li, Lj, sigFigs);
            Vi(ii,jj) = iMonge;
            iMonge = iMonge + 1;
        end
    end

    % Symmetrize the index matrix        
    Vi = Vi + Vi';
end

function V = makeMonge(N, M, sigFigs)
    % Start by constructing a square matrix for simplicity
    maxL = max(N, M);    
    while true
        % Monotonic sequence
        %     Vm(:,:,2) = Vm(:,:,1).^2; % L2 term
        seq1 = cumsum(rand(1, maxL));
        seq2 = cumsum(rand(1, maxL));
        V = abs(bsxfun(@minus, seq1, seq2'));
        if sigFigs
            V = fround(V, sigFigs);
        end
        if IsMonge(V)            
            break;
        end        
    end
    
    % Slice to fit
    V = V(1:N,1:M);
    
    assert(IsMonge(V));
end

%-------------------------------------------------------------------------%
function tf = IsMonge(M)
%
% Checks if matrix M is a Monge matrix.
%
% The following must hold for a Monge matrix (with finite entries):
%   M_ik + M_jl <= M_il + M_jk \forall 1 <= i <= j <= m, 1 <= k <= l <= n
%
% M may be 3D array, in which case each 2D "slice" is checked for Monge
% property
%

tf = true;

for si=1:size(M,3)
    [m n] = size(M(:,:,si));
    for ii=1:(m-1)
        for kk=1:(n-1)
            tf =  ( M(ii,kk) + M(ii+1,kk+1) <= M(ii, kk+1) + M(ii+1, kk) );
            if ~tf, return; end
        end
    end
end

end
