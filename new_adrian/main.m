%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);

%% Rest of it
N = 4;
%edgeProb=0.5;
edgeProb=1;

seed = -1;

%epsilon = 0.05;
%epsilon = 0.01;
%epsilon = 0.1;

Tmax=0;  % Local potentials
%Wmax=5;

%assoc=0;
assoc=1;

maxiter=100;

methods = {'simple', 'minsum', 'adaptivesimple', 'adaptiveminsum'};

nTrials = 10;

clear N1Hist H2Hist N2_2Hist fdHist;
N1Hist(nTrials) = 0;
N2Hist(nTrials) = 0;
N2_2Hist(nTrials) = 0;
fdHist(nTrials, 3) = 0;

maxPoints = 1e4;

for t = 1:nTrials
    [ gamma1,gamma2,gamma2_2,N1,N2,N2_2,seed,Am,Bm,theta,W,K,zeta, J, thisN1,thisN2,thisN2_2,L,U ] = gamma12( N,epsilon,Tmax,Wmax,edgeProb,assoc,seed,maxiter );

    if any(1 - Bm - Am == 0)
        fprintf(1, 'Problem is too trivial; moving on');
        continue;
    end
    
    W = sparse(W);
    
    [lbpLogZ, ~, lbpOneMarg, lbpTwoMarg, ~] = solveDAI(theta, W, 'BP', '[inference=SUMPROD,updates=SEQFIX,logdomain=0,tol=1e-9,maxiter=10000,damping=0.0]');
    [logZ1True, ~, ~, ~, ~] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');        

    betheOpts = struct('useMooij', true, 'bbThresh', 1e-4, 'maxIter', 200);
    [oldBethe, bOneMarg, bTwoMarg, bMisc] = BetheApprox_opt_mex(theta, W, epsilon, betheOpts);

    fprintf(1, 'True logZ = %g ; LBP logZ = %g ; Bethe logZ = %g\n', logZ1True, lbpLogZ, oldBethe);
%    logZ1True = 0;

    %% One run
    if N1 < maxPoints
        gams1 = fillGams(gamma1, Am, Bm, true);
        
        %[logZ1, oneMarginals, twoMarginals, misc] = solveBetheNew(theta, W, gams1);
        [logZ1, oneMarginals, twoMarginals, misc] = BetheGams_mex(theta, W, gams1);        
        logZ1Gap = abs(logZ1 - lbpLogZ);
        fprintf(1, 'LOGZ1 GAP: %g\n', abs(logZ1Gap));
        if logZ1Gap > epsilon
            warning('gamma1 bounds failed');
        end
        
        N1Hist(t) = N1;
    else
        warning('Method 1 skipped because of %d points', N1);
        N1Hist(t) = -1;
    end
    
    if N2 < maxPoints
        % Note that we don't actually run this code because it's the same as the old Bethe code.
        % Although we count the number of generated points anyways.
%        gams2 = fillGams(gamma2, Am, Bm, true);
%        
%        [logZ2, oneMarginals, twoMarginals, misc] = BetheGams_mex(theta, W, gams2);
%        logZ2Gap = abs(logZ2 - lbpLogZ);
%        fprintf(1, 'LOGZ2 GAP: %g\n', logZ2Gap);
%        if logZ2Gap > epsilon
%            warning('gamma2 bounds failed');
%        end
        
        N2Hist(t) = N2;
    else
        warning('Method 2 skipped because of %d points', N2);
        N2Hist(t) = -1;
    end
    
    if N2_2 < maxPoints
        gams2_2 = fillGams(gamma2_2, Am, Bm, true);    
        
        [logZ2_2, oneMarginals, twoMarginals, misc] = BetheGams_mex(theta, W, gams2_2);        
        logZ2_2Gap = abs(logZ2_2 - lbpLogZ);
        fprintf(1, 'LOGZ2_2 GAP: %g\n', logZ2_2Gap); 
        if logZ2_2Gap > epsilon
            warning('gamma2_2 bounds failed');
        end
        
        N2_2Hist(t) = N2_2;
    else
        warning('Method 2_2 skipped because of %d points', N2_2);
        N2_2Hist(t) = -1;
    end    
    
    gams = cell(N, 3);
    logZ_fd = zeros(3,1);
    
    for im = 1:4 
        method = methods{im};
        [ gm,nPts ] = fdm(theta, W, Am, Bm, epsilon, method, L, U);

        if nPts > maxPoints
            warning('First derivative method %d skipped because of %d points', im, nPts);
            continue;
        end
            
        if im == 1 || im == 2
            assert(all(~isnan(gm)));
            gams(:,im) = fillGams(gm, Am, Bm, false);
        else
            assert(iscell(gm));
            gams(:,im) = gm;
        end

        [logZ2_fd(im), oneMarginals, twoMarginals, misc] = BetheGams_mex(theta, W, gams(:,im));
        gap = abs(logZ2_fd(im) - lbpLogZ);
        fprintf(1, 'LOGZ2_fd GAP: %g\n', gap);
        if gap > epsilon
            warning('adaptive method %d failed\n');
        end
        
        fdHist(t,im) = nPts;
    end
end
    

