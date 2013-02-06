function showElMat( elMat )

    maxI = max(elMat(:,1));
    maxJ = max(elMat(:,2));
    nNodes = max(maxI, maxJ) - 2;    
    
    sNode = nNodes + 1;
    tNode = nNodes + 2;
    
    iVec = elMat(:,1);
    jVec = elMat(:,2);
    ijVec = elMat(:,3);
    
    connMat = sparse(iVec, jVec, ijVec, nNodes + 2, nNodes + 2);
    
    gObj = biograph(connMat, []);

    figure;
    view(gObj);    

end

