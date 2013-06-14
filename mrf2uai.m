function mrf2uai(mrf, outFile)
% mrf2uai  Write mrf to UAI-formatted file.
%
%   See http://www.cs.huji.ac.il/project/PASCAL/fileFormat.php
%   mrf format is a structure with fields are arguments to 
%   MultiLabelSubModular

    % All factors here are pairwise!
    [siVec, sjVec, swVec] = findUT(mrf.newW);    
    nNodes = length(mrf.D);
    nEdges = length(siVec);
    
    % Convert the energy functional to a factor graph.
    fd = fopen(outFile, 'w');
    fprintf(fd, 'MARKOV\n');                         % 1: must say MARKOV
    fprintf(fd, '%d\n', nNodes);                     % 2: number of nodes
    fprintf(fd, '%d ', cellfun(@length, mrf.D));     % 3: cardinalities
    fprintf(fd, '\n%d\n', nNodes + nEdges);          % 4: number of factors
    
    % Each line offset + i specifies the SCOPE of factor i. First entry is
    % the number of nodes in the factor; rest are labels of the factor.    
    %
    % Note that UAI uses ZERO indices.        
    ab = [ones(1, nNodes); (1:nNodes) - 1];
    fprintf(fd, '%d %d\n', ab(:));
    
    % Interleave siVec and sjVec. They are currently column vectors. This
    % results in printing a list of edges.    
    edges = [2 * ones(1, nEdges) ; siVec' - 1 ; sjVec' - 1];
    fprintf(fd, '%d %d\n', edges(:));
    
    % Now, print the factor tables. The first line is the size of the
    % table, and then we print the table.
    %
    % We stored the energies, but MPLP expects log probabilities. So
    % negate.
    for n = 1:nNodes
        fprintf(fd, '%d\n', length(mrf.D{n}));
        fprintf(fd, '%d ', -mrf.D{n});        
        fprintf(fd, '\n');
    end
    
    fmt = get(0, 'Format');
    format long;
    for ne = 1:nEdges
        fprintf(fd, '%d\n', numel(mrf.Vm{ne}));
        fprintf(fd, '%s\n', toString(mrf.Vm{ne}, 'disp'));
    end
    format(fmt);
    
    fclose(fd);
end

