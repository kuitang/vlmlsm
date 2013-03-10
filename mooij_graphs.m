%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);

N = 4;
j = 1;
t = 0;

J = j * ones(N,N);
J = J - diag(diag(J));
eta = t * ones(N,1);

epsilon = 1e-1;

% convert from (-1,+1) to (0,1) convention
theta = 2 * eta - 2 * sum(J,2); % zeta is called theta in the paper
W = sparse(4 * J);

opts = struct();

[trueLogZ, ~, trueOneMarg, trueTwoMarg, JTTime] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');           
disp(['JTree Finished']);
[lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, lbpTime] = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');
disp(['LBP Finished']);

fprintf(1, 'logZ: True: %g; LBP: %g\n', trueLogZ, lbpLogZ);

% [logZ, oneMarg, twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, opts);   
% disp(['Bethe Finished']);

[A, B, ~] = BBP(theta, W, 0.002, 1000);
[isz, lambda, theorybound] = getIntervalSz(A, B, W, 1e-2)
1 - B - A

[Am, Bm] = bpbound(N, theta, W, 1000);
[iszm, lambda, theoryboundm] = getIntervalSz(Am, Bm, W, 1e-2)
1 - Bm - Am

%'BP[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]'

%fprintf(1, 'logZ: True: %g, Bethe: %g ; LBP: %g\n', trueLogZ, logZ, lbpLogZ);
