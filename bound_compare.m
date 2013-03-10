%% Setup
path_to_dai = '../libDAI-0.3.1/matlab';
addpath(path_to_dai);
nTrials = 100;

%% Problem
%nNodes = 10000;
nNodes = 16;
%testTrees = true;
testTrees = false;
treeLoops = 200;

boundIters = 200;

A(nNodes,nTrials) = 0;
B(nNodes,nTrials) = 0;
ranges(nNodes,nTrials) = 0;

Am(nNodes,nTrials) = 0;
Bm(nNodes,nTrials) = 0;
mRanges(nNodes,nTrials) = 0;

% Uniform bounds for weights
ta = -2;
tb = 2;
wb = 1;

density = 1;

epsilon = 1e-2;

for t = 1:nTrials
    %% Set up the problem
    if testTrees
        T = randTree(nNodes, nLoops);
        [theta, W] = makeUnifProblem(nNodes, T, ta, tb, wb )
    else
        [theta, W] = makeUnifProblem(nNodes, density, ta, tb, wb);
    end
    
    % Mooij's bounds
    [Am(:,t), Bm(:,t)]  = bpbound(nNodes, theta, W, boundIters);
    [A(:,t), B(:,t), ~] = BBP(theta, W, 0, boundIters);
    
    mRanges(:,t) = 1 - Bm(:,t) - Am(:,t);
    ranges(:,t)  = 1 - B(:,t) - A(:,t);
    
    fprintf(1, 'Mean and max interval length: Mooij: %g %g Mine: %g %g\n', ...
            mean(mRanges(:,t)), max(mRanges(:,t)), mean(ranges(:,t)), max(ranges(:,t)));    
end

figure;
hist(mean(mRanges, 1), 100);
title('Mean Mooij interval lengths');

figure;
hist(max(mRanges, 1), 100);
title('Max Mooij interval lengths');

figure;
hist(mean(ranges, 1), 100);
title('Mean BBP interval lengths');

figure;
hist(max(ranges, 1), 100);
title('Max BBP interval lengths');
