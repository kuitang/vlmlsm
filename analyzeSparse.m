%% Setup part:

nNodeVec = [10 30 50];
nFiles   = length(nNodeVec);

commaStr = [];
for n = 1:nFiles
    nodeLabels{n} = num2str(nNodeVec(n));
    blankLabels{n} = '';
    commaStr = [commaStr ', ' nodeLabels{n}];
end

commaStr(1:2) = [];

% Assume uniform across files
nTrials = 100;

%% Collate results from .mat files
kappas(nTrials,nFiles) = 0;
theoryBounds(nTrials,nFiles) = 0;

allBetheTimes(nTrials,nFiles) = 0;
allMooijTimes(nTrials,nFiles) = 0;
allSumBetheGap(nTrials,nFiles) = 0;
allSumMooijGap(nTrials,nFiles) = 0;

totDiffs = cell(nFiles, 1);

%%
for n = 1:nFiles
    nNodes = nNodeVec(n);
    fn     = ['sparse_results/' num2str(nNodes) '_e01.mat'];
    load(fn);
    
    % Collate the data we want
    
    % For the "Average error of one-marginals"
    totDiffs{n} = totDiff;
    
    % Prepare our actual runtimes vs the theoretical bounds
    for t = 1:nTrials
        [~, kappas(t,n), theoryBounds(t,n)] = getIntervalSz(A(:,t),B(:,t), W, epsilon);
    end
    
    allBetheTimes(:,n) = betheTimes(:,1);
    allMooijTimes(:,n) = betheTimes(:,2);
    
    allSumBetheGap(:,n) = sum(ABgap(:,:,1),1)';
    allSumMooijGap(:,n) = sum(ABgap(:,:,2),1)';
end

save sparse_collated.mat allBetheTimes allMooijTimes allSumBetheGap allSumMooijGap kappas theoryBounds;

%% Plot errors -- IN PUBLICATION AS OF 13 MARCH
epsilonTxt = [' for \epsilon = ' num2str(epsilon)];

fig('units','inches','width',6.5,'height',1.75,'font','Helvetica','fontsize',10);
for n = 1:nFiles
    subplot(1, nFiles, n);    
    hist(totDiffs{n}(:,1));
    title(['Avg. err, n = ' num2str(nNodeVec(n))]);
end

%% Plot runtime vs nodes -- IN PUBLICATON AS OF 13 MARCH
fig('units','inches','width',3,'height',3,'font','Helvetica','fontsize',10);
hold on;
% TODO: Remove hardcoded labels
boxplot(allBetheTimes, 'labels', nodeLabels, 'colors', 'g');
boxplot(theoryBounds,  'labels', blankLabels, 'colors', 'b');
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


%% Plot runtime vs Kappa -- IN PUBLICATION AS OF 13 MARCH

fig('units','inches','width',3,'height',3,'font','Helvetica','fontsize',10);
hold on;
colors = {'r', 'm', 'b'};

for n = 1:nFiles        
    [xs, ix] = sort(kappas(:,n));
    ys = allBetheTimes(:,n);
    ys = ys(ix);
    ts = theoryBounds(:,n);
    ts = ts(ix);

    scatter(xs, ys, colors{n}, 'x');
    set(gca, 'XScale', 'Log');
    set(gca, 'YScale', 'Log');
    
    plot(xs, ts, 'color', colors{n});
    xlabel('\kappa');
    ylabel('Runtime');    
end

title(['Runtime vs \kappa; n = ' commaStr]);

hold off;

%% Plot BBP Width / Mooij Width for the 16 node problem
% Implicitly rely on the last instance saved in ABgap
gapRatio = ABgap(:,:,1) ./ ABgap(:,:,2);

%% Mooij and Bethe vs Kappa
figure;
[xs, ix] = sort(kappas(:,end));
bs = allBetheTimes(:,end);
ms = allMooijTimes(:,end);
hold on;
scatter(xs, ms, 'g', 'x');
scatter(xs, bs, 'r', 'x');
set(gca, 'XScale', 'Log');
set(gca, 'YScale', 'Log');
legend('Mooij', 'Bethe');
title(['Times vs. \kappa, n = ' num2str(nNodeVec(end))]);

% NOTE: By inspection, intervalsize is NOT interesting; it varies by only
% half its magnitude. (Thus, the curvature (Kappa) bound seems to account 
% for both a slightly smaller interval and more intervals?)

%% BBP vs. Mooij widths
gapRatios(nTrials) = 0;
for t = 1:nTrials
    gapRatios(t) = sum(vec(ABgap(:,t,1))) / sum(vec(ABgap(:,t,2)));
end

figure
hist(gapRatios);
title(['AB gap ratio (BBP / Mooij) over ' num2str(nTrials) ' trials']);

% Time tradeoff plot
figure;
[xs, ix] = sort(betheTimes(:,1));
plot(xs, gapRatios(ix));
title('BBP/Mooij width vs. Bethe time');
set(gca, 'XScale', 'Log');
