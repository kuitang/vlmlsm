%% Setup part:

nNodeVec = [4 8 16];
nFiles   = length(nNodeVec);

commaStr = [];
for n = 1:nFiles
    nodeLabels{n} = num2str(nNodeVec(n));
    commaStr = [commaStr ', ' nodeLabels{n}];
end

commaStr(1:2) = []

% Assume uniform across files
nTrials = 100;

%% Collate results from .mat files
lambdas(nTrials,nFiles) = 0;
theoryBounds(nTrials,nFiles) = 0;

allBetheTimes(nTrials,nFiles) = 0;
allMooijTimes(nTrials,nFiles) = 0;
allSumBetheGap(nTrials,nFiles) = 0;
allSumMooijGap(nTrials,nFiles) = 0;

totDiffs = cell(nFiles, 1);

%%
for n = 1:nFiles
    nNodes = nNodeVec(n);
    fn     = [num2str(nNodes) '.mat'];
    load(fn);
    
    % Collate the data we want
    
    % For the "Average error of one-marginals"
    totDiffs{n} = totDiff;
    
    % Prepare our actual runtimes vs the theoretical bounds
    for t = 1:nTrials
        [~, lambdas(t,n), theoryBounds(t,n)] = getIntervalSz(A(:,t),B(:,t), W, epsilon);
    end
    
    allBetheTimes(:,n) = betheTimes(:,1);
    allMooijTimes(:,n) = betheTimes(:,2);
    
    allSumBetheGap(:,n) = sum(ABgap(:,:,1),1)';
    allSumMooijGap(:,n) = sum(ABgap(:,:,2),1)';
end

save collated.mat allBetheTimes allMooijTimes allSumBetheGap allSumMooijGap lambdas theoryBounds;

%% Plot errors
epsilonTxt = [' for \epsilon = ' num2str(epsilon)];

fig('units','inches','width',6.5,'height',1.75,'font','Helvetica','fontsize',10);
for n = 1:nFiles
    subplot(1, nFiles, n);    
    hist(totDiffs{n}(:,1));
    title(['Avg. err, n = ' num2str(nNodeVec(n))]);
end

%% Plot runtime vs nodes
fig('units','inches','width',3,'height',3,'font','Helvetica','fontsize',10);
hold on;
% TODO: Remove hardcoded labels
boxplot(allBetheTimes, 'notch', 'marker', 'labels', {'4', '8', '16'}, 'colors', 'g');
boxplot(theoryBounds, 'notch', 'marker',  'labels', {'', '', ''}, 'colors', 'b');
set(gca, 'YScale', 'Log');
title(['Runtime and bound (top) vs nodes; ' num2str(nTrials) ' trials']);
hold off;

% fig('units','inches','width',3,'height',3,'font','Helvetica','fontsize',12);
% hold on;
% boxplot(allMooijTimes, 'notch', 'marker', 'labels', nodeLabels, 'colors', 'g');
% boxplot(theoryBounds, 'notch', 'marker',  'labels', nodeLabels, 'colors', 'b');
% set(gca, 'YScale', 'Log');
% title(['Mooij runtime (bottom) and theoretical bound (top) vs nodes; ' num2str(nTrials) ' trials']);
% hold off


%% Plot runtime vs Lambda

fig('units','inches','width',3,'height',3,'font','Helvetica','fontsize',10);
hold on;
colors = {'r', 'm', 'b'};

for n = 1:nFiles        
    [xs, ix] = sort(lambdas(:,n));
    ys = allBetheTimes(:,n);
    ys = ys(ix);
    ts = theoryBounds(:,n);
    ts = ts(ix);

    scatter(xs, ys, colors{n}, 'x');
    set(gca, 'XScale', 'Log');
    set(gca, 'YScale', 'Log');
    
    plot(xs, ts, 'color', colors{n});
    xlabel('Lambda');
    ylabel('Runtime');    
end

title(['Runtime vs Lambda; n = ' commaStr]);

hold off;
% NOTE: By inspection, intervalsize is NOT interesting; it varies by only
% half its magnitude. (Thus, the curvature (Lambda) bound seems to account 
% for both a slightly smaller interval and more intervals?)
