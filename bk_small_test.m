%% CLRS 726
disp('Running example on CLRS pp. 726 (3rd ed)');
nNodes = 4;
s = -1;
t = -2;
iVec  = [1 2 3 2 4 s s 3 4]';
jVec  = [3 1 2 4 3 1 2 t t]';
ijVec = [12 4 9 14 7 16 13 20 4]';
jiVec = zeros(length(ijVec), 1);

[e, cut] = BK_mex(iVec, jVec, ijVec, jiVec, nNodes, s, t)
% Looking at the picture, the cut with v3 and t touches the -> t edges
% 12, 7, 4, adding up to 23.
assert(all(cut(2:end) == [0 0 1 0]));
assert(e == 23);

%% CLRS 728
disp('Running example on CLRS pp. 728 (3rd ed)');
nNodes = 2;
iVec = [s s 1 1 2];
jVec = [1 2 2 t t];
ijVec = [1e6 1e6 1 1e6 1e6];
jiVec = zeros(1, length(ijVec));

[e, cut] = BK_mex(iVec, jVec, ijVec, jiVec, nNodes, s, t)
assert(e == 2e6);

