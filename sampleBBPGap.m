%% Parameters
nNodes = 16;
nTrials = 1000;
epsilon = 1e-2;

A = zeros(nNodes, nTrials);
B = zeros(nNodes, nTrials);
theta = zeros(nNodes, nTrials);
gap = zeros(nNodes, nTrials);
W = cell(1, nTrials);

%% Run

for t = 1:nTrials
    [theta(:,t), W{t}] = makeProblem(nNodes);
    [A(:,t), B(:,t), ~] = BBP(theta(:,t), W{t}, 0.002, 200);
    intervalSz = getIntervalSz(A(:,t), B(:,t), W{t}, epsilon);
    gap(:,t) = 1 - B(:,t) - A(:,t);
    
    assert(all(gap(:,t) > 0), 'Gap must be positive.');
    
    totGap(t) = sum(gap(:,t));
    stdGap(t) = std(gap(:,t));
    intervalSzs(t) = intervalSz;
    nIntervals(t) = floor(totGap(t) / intervalSz) + 2 * nNodes;
end

figure;
hist(totGap, 100);
title('Total gap size');

figure;
hist(stdGap, 100);
title('Gap standard deviations');

figure;
hist(intervalSzs, 100);
title('Interval sizes');

figure;
hist(nIntervals, 100);
title('Total number of intervals');


