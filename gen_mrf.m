function [ D, W, Vi, Vm ] = gen_mrf( N, maxL, edge_prob, figs )
% gen_mrf Generate test cases for MultiLabelSubModular
%   [D, W, Vi, Vm] = gen_mrf(N, L, Nmonge, edge_prob)
%
%   N         - number of nodes
%   maxL      - Maximum number of labels
%   edge_prob - of all N*(N-1)/2 edges, what fraction will appear?
%   figs      - if present, round to figs places after the decimal
%               (generate nice numbers for debug cases)

    if nargin <= 3
        figs = 0; % Sentinel
    end
    
    D = cell(N, 1);
    for n = 1:N
        % Must have at least two states!
        L = 1 + randi(maxL - 1, 1);
        D{n} = rand(1, L);
        if figs
            D{n} = fround(D{n}, figs);
        end
    end    
    
    nzmat = 0;
    while sum(nzmat(:)) == 0
        nzmat = rand(N, N) < edge_prob;
        % remove diagonal;       
        % http://www.mathworks.com/matlabcentral/newsreader/view_thread/47182
        nzmat(1:N+1:N*N) = 0;   
    end
    
    %W = ones(N, N) .* nzmat;        
    W = abs(rand(N, N) .* nzmat);
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
            
            Vm{iMonge} = makeMonge(Li, Lj, figs);
            Vi(ii,jj) = iMonge;
            iMonge = iMonge + 1;
        end
    end

    % Symmetrize the index matrix        
    Vi = Vi + Vi';
end
