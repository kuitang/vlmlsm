% Make instances of hard problems (the interior of the top triangle in
% Mooij and Kappen 2005).
%
% LBP *will* fail in this regime, but both BBP and Mooij-Kappen bounds will
% be very loose, and interval sizes for both problems will be tiny. Expect
% long runtime.

%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 100;

nNodes = 4;
        
printStats = true;

% trueLogZ
epsilon = 0.1;

oneMarg(nNodes,nTrials) = 0;
trueOneMarg(nNodes,nTrials) = 0;
lbpOneMarg(nNodes,nTrials) = 0;
A(nNodes,nTrials) = 0;
B(nNodes,nTrials) = 0;
ABgap(nNodes,nTrials) = 0;
intervalSzs(nTrials) = 0;

% TODO: Test various sparsity patterns.

%% Zero-out some vectors... does it even matter?
totDiff(nTrials) = 0;
lbpTotDiff(nTrials) = 0;

trueLogZGaps(nTrials) = 0;
lbpLogZGaps(nTrials) = 0;

betheTimes(nTrials) = 0;
bkEdges(nTrials) = 0;
maxFlowTimes(nTrials) = 0;

lbpTimes(nTrials) = 0;
JTTimes(nTrials) = 0;

mooijOpts = struct('useMooij', true);

problems = cell(1, nTrials);

successVec(nTrials) = false;

fn = ['results ' num2str(nNodes) ' nodes ' datestr(now, 0);];
%% Loop
for t = 1:nTrials
    
%     % Draw (t,j) from the triangle (-1,2), (1,2), (0,1)
%     j  = unifrnd(1, 2)                          % Uniform coupling strength            
%     dt = j - 1;
%     th = unifrnd(-dt, dt)
%     
%     J = j * ones(nNodes,nNodes);
%     J = J - diag(diag(J));
%     eta = th * ones(nNodes,1);  

    % Random graphs, instead of uniform potentials. Zero local potential.
    j = 0;
    eta = j * ones(nNodes, 1);
    sigma = 0.9;
    
    J = abs(normrnd(0.5 + j, sigma, nNodes, nNodes));
    % Symmetrize;
    J = 0.5 * (J + J');
    % Remove diagonal
    J = J - diag(diag(J));

    % convert from (-1,+1) to (0,1) convention
    theta = 2 * eta - 2 * sum(J,2); % zeta is called theta in the paper
    W = sparse(4 * J);

    problems{t} = struct('eta', eta, 'theta', theta, 'W', W);

    if printStats
        [A, B, ~] = BBP(theta, W, 0.002, 1000)
        [isz, kappa, theorybound] = getIntervalSz(A, B, W, 1e-2);
        1 - B - A;
        gaps = 1 - B - A;
        totalGap = sum(gaps(:));
        nIntervals = totalGap / isz

        [Am, Bm] = bpbound(nNodes, theta, W, 1000)
        [iszm, kappa, theoryboundm] = getIntervalSz(Am, Bm, W, 1e-2)
        mGaps = 1 - Bm - Am;
        totalMGaps = sum(mGaps(:));
        nMIntervals = totalMGaps / iszm

        nIntervalRatio = nMIntervals / nIntervals
    end
    
    %% Run JTree and LBP
    [trueLogZ, ~, trueOneMarg(:,t), trueTwoMarg, JTTimes(t)] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');           
    disp(['JTree Finished']);
    [lbpLogZ, ~, lbpOneMarg(:,t), lbpTwoMarg, lbpTime(t)] = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');
    disp(['LBP Finished']);
    fprintf(1, 'logZ: True: %g; LBP: %g\n', trueLogZ, lbpLogZ);       
    
    
    %% Run the time-killer
    try
      [mooijlogZ, oneMarg(:,t,2), twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, mooijOpts);            
      fprintf(1, 'logZ: Mooij: %g; True: %g; LBP: %g\n', mooijlogZ, trueLogZ, lbpLogZ);       
      successVec(t) = true;
    catch err
      disp(err)
      successVec(t) = false;
    end
    
    % [logZ, oneMarg, twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, opts);   
    % disp(['Bethe Finished']);



end
