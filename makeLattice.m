function [ W ] = makeLattice( rows, cols )
% makeLattice Make logical 4-neighbor 2d lattice
%
% W = makeLattice(rows, cols, eight) makes a 2d lattice of rows x cols.
% If eight is nonempty, makes 8-neighbor lattice. Else, four-neighbor.
%
% W can be used as sparsity pattern for randomly testing MRF codes.

    dims = [rows cols];        
    
    nodeMax = 8 + (cols - 2) * (rows - 2);
    iVec(nodeMax) = 0;
    jVec(nodeMax) = 0;
    
    % If we use 4-neighbors, we need neighbors explicitly.
    iVec(1:2) = sub2ind(dims, 1, 1);
    jVec(1)   = sub2ind(dims, 1, 2);
    jVec(2)   = sub2ind(dims, 2, 1);

    iVec(3:4) = sub2ind(dims, 1, cols);
    jVec(3)   = sub2ind(dims, 1, cols - 1);
    jVec(4)   = sub2ind(dims, 2, cols);

    iVec(5:6) = sub2ind(dims, rows, 1);
    jVec(5)   = sub2ind(dims, rows - 1, 1);
    jVec(6)   = sub2ind(dims, rows, 2);

    iVec(7:8) = sub2ind(dims, rows, cols);
    jVec(7)   = sub2ind(dims, rows - 1, cols);
    jVec(8)   = sub2ind(dims, rows, cols - 1);

    ind = 9;
    
    for j = 2:(cols-1)
        for i = 2:(rows-1)
            iVec(ind:(ind+3)) = sub2ind(dims, i, j);
            jVec(ind)         = sub2ind(dims, i, j+1);
            jVec(ind+1)       = sub2ind(dims, i+1, j);
            jVec(ind+2)       = sub2ind(dims, i-1, j);
            jVec(ind+3)       = sub2ind(dims, i, j-1);
            ind = ind + 4;            
        end
    end
    
    wVec = true(size(iVec));
    W = sparse(iVec, jVec, wVec);    
end

