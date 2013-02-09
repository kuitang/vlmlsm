function [ siVec, sjVec, swVec ] = findUT( W )
% Find the upper triangular of a (sparse) W.

    [iVec, jVec, wVec] = find(W);
    % As always, take the upper triangular
    sel = iVec < jVec;
    siVec = iVec(sel);
    sjVec = jVec(sel);
    swVec = wVec(sel);

end

