%% Heskes pp. 2408
%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nNodes = 9;
% Alpha is critical at 2/3, ~0.79, and ~0.88 (LBP will diverge for higher)
% NOTE: Subract from one here.
alpha = 0.1;
% Reparameterization to our (Eq 1), up to constant alpha
theta = (1 - 2*alpha) * ones(nNodes, 1);
w = 4*alpha - 2;

epsilon = 1e-3;
opts = struct();

% Toroid
edges = [...
1 2
1 3
1 4
1 7
2 3
2 5
2 8
3 6
3 9
4 5
4 7
4 6
5 6
5 8
6 9
7 8
7 9
8 9 ];

W = sparse(vertcat(edges(:,1), edges(:,2)), ...
           vertcat(edges(:,2), edges(:,1)), ...
           w * ones(2 * size(edges, 1), 1), ...
           nNodes, nNodes);

%% blah blah
[betheLogZ, betheOneMarg, betheTwoMarge, betheMisc] = BetheApprox_opt_mex(theta, W, epsilon, opts);
[trueLogZ, trueJoint, trueOneMarg, trueTwoMarg, JTTime]  = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');
[lbpLogZ,  lbpJoint,  lbpOneMarg,  lbpTwoMarg,  lbpTime] = solveDAI(theta, W, 'BP', '[tol=1e-9,logdomain=1,updates=SEQRND]');

[betheOneMarg lbpOneMarg trueOneMarg]
%[lbpOneMarg trueOneMarg]

% Lattice
edges = makeLattice(3, 3);

%% Lattice with two removed edges
edges = makeLattice(3, 3);
edges(1,4) = 0;
edges(4,1) = 0;
edges(6,9) = 0;
edges(9,6) = 0;
 
