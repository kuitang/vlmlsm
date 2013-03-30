% Make instances of hard problems (the interior of the top triangle in
% Mooij and Kappen 2005).
%
% LBP *will* fail in this regime, but both BBP and Mooij-Kappen bounds will
% be very loose, and interval sizes for both problems will be tiny. Expect
% long runtime.

% nTrials and randGraph set from the environment!

%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nNodes = 5;
        
printStats = true;
runBethe = true;
runMooij = false;

% trueLogZ
epsilon = 10;

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

% 1 - true, 2 - lbp, 3 - mooij, 4 - bethe
logZ(nTrials, 4) = 0;

betheTimes(nTrials) = 0;
bkEdges(nTrials) = 0;
maxFlowTimes(nTrials) = 0;

lbpTimes(nTrials) = 0;
JTTimes(nTrials) = 0;

mooijOpts = struct('useMooij', true);
betheOpts = struct('useMooij', false);

nIntervalCap = 1e5;

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
    if randGraph
        % Generate random graph according to MK Fig 5, but flipped to
        % positive edge weights
        accepted = false;
        while ~accepted
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


            if printStats            
                %% Bethe stats
                [A, B, ~] = BBP(theta, W, 0.002, 1000)
                [isz, kappa, theorybound] = getIntervalSz(A, B, W, epsilon);
                1 - B - A;
                gaps = 1 - B - A;
                totalGap = sum(gaps(:));
                nIntervals = totalGap / isz;
                
                %% Mooij stats
                [Am, Bm] = bpbound(nNodes, theta, W, 1000)
                [iszm, kappa, theoryboundm] = getIntervalSz(Am, Bm, W, epsilon)
                mGaps = 1 - Bm - Am;
                totalMGaps = sum(mGaps(:));
                nMIntervals = totalMGaps / iszm;

                nIntervalRatio = nMIntervals / nIntervals;
            end

            if nMIntervals < nIntervalCap && nIntervals < nIntervalCap
                accepted = true;
                disp('Interval gap reached. Proceeding.');
            end
        end
    else
        %% Use the (eta, j) = (0, 1) graph
        eta = 0;
        J   = ones(nNodes);
        J   = J - diag(diag(J));
        
        % convert from (-1,+1) to (0,1) convention
        theta = 2 * eta - 2 * sum(J,2); % zeta is called theta in the paper
        W = sparse(4 * J);        
    end
    
    problems{t} = struct('eta', eta, 'theta', theta, 'W', W);
    
    %% Run JTree and LBP
    [logZ(t,1), ~, trueOneMarg(:,t), trueTwoMarg, JTTimes(t)] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');           
    disp(['JTree Finished']);
    [logZ(t,2), ~, lbpOneMarg(:,t), lbpTwoMarg, lbpTime(t)] = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=PARALL,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');
    disp(['LBP Finished']);
    %fprintf(1, 'logZ: True: %g; LBP: %g\n', logZ(t,1), logZ(t,2));       
    
    %% Run the time-killer
    try
      if runMooij
        [logZ(t,3), oneMarg(:,t,2), twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, mooijOpts);            
      end
      if runBethe
        [logZ(t,4), oneMarg(:,t,2), twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, betheOpts);
      end
      fprintf(1, 'logZ: True: %g, LBP: %g, Mooij: %g; Bethe: %g\n', logZ(t,:));
      successVec(t) = true;
    catch err
      disp(err)
      successVec(t) = false;
    end
    
    % [logZ, oneMarg, twoMarg, misc] = BetheApprox_opt_mex(theta, W, epsilon, opts);   
    % disp(['Bethe Finished']);
end

save(fn);

