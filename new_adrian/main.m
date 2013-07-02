%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);

%% Rest of it
N = 4;
%edgeProb=0.5;
edgeProb=1;

seed = -1;

%epsilon = 0.05;
epsilon = 0.01;

Tmax=0;  % Local potentials
Wmax=10;

%assoc=0;
assoc=1;

maxiter=100;

methods = {'simple', 'lagrangian', 'adaptive' };

nTrials = 10;

clear N1Hist H2Hist N2_2Hist fdHist;
N1Hist(nTrials) = 0;
N2Hist(nTrials) = 0;
N2_2Hist(nTrials) = 0;
fdHist(nTrials, 3) = 0;

maxPoints = 1e4;

for t = 1:nTrials
    [ gamma1,gamma2,gamma2_2,N1,N2,N2_2,seed,Am,Bm,theta,W,K,zeta, J, thisN1,thisN2,thisN2_2,L,U ] = gamma12( N,epsilon,Tmax,Wmax,edgeProb,assoc,seed,maxiter );
    
%     [logZ1True, ~, ~, ~, ~] = solveDAI(theta, W, 'JTREE', '[updates=HUGIN,verbose=0]');        
    logZ1True = 0;

    %% One run
    if N1 < maxPoints
        gams1 = fillGams(gamma1, Am, Bm, true);
        
        [logZ1, oneMarginals, twoMarginals, misc] = solveBetheNew(theta, W, gams1);
        fprintf(1, 'LOGZ1 GAP: %g\n', abs(logZ1 - logZ1True));
        
        N1Hist(t) = N1;
    else
        warning('Method 1 skipped because of %d points', N1);
        N1Hist(t) = -1;
    end
    
    if N2 < maxPoints
        gams2 = fillGams(gamma2, Am, Bm, true);
        
        [logZ2, oneMarginals, twoMarginals, misc] = solveBetheNew(theta, W, gams2);
        fprintf(1, 'LOGZ2 GAP: %g\n', abs(logZ2 - logZ1True));
        
        N2Hist(t) = N2;
    else
        warning('Method 2 skipped because of %d points', N2);
        N2Hist(t) = -1;
    end
    
    if N2_2 < maxPoints
        gams2_2 = fillGams(gamma2_2, Am, Bm, true);    
        
        [logZ2_2, oneMarginals, twoMarginals, misc] = solveBetheNew(theta, W, gams2_2);        
        fprintf(1, 'LOGZ2_2 GAP: %g\n', abs(logZ2_2 - logZ1True));
        
        N2_2Hist(t) = N2_2;
    else
        warning('Method 2_2 skipped because of %d points', N2_2);
        N2_2Hist(t) = -1;
    end    
    
    gams = cell(N, 3);
    logZ_fd = zeros(3,1);
    
    for im = [1 3]
%    for im = 1:3
        method = methods{im};
        [ gm,nPts ] = fdm(theta, W, Am, Bm, epsilon, method, L, U);

        if im ~= 3            
            assert(all(~isnan(gm)));
            if nPts < maxPoints
                warning('First derivative method %d skipped because of %d points', im, nPts);
            else
                im
                gams(:,im) = fillGams(gm, Am, Bm, false);

                [logZ2_fd(im), oneMarginals, twoMarginals, misc] = solveBetheNew(theta, W, gams(:,im));
                fprintf(1, 'LOGZ2_fd GAP: %g\n', abs(logZ2_fd(im) - logZ1True));
            end
        end
        
        fdHist(t,im) = nPts;
    end

end
    

