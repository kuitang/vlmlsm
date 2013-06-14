%% What size for the publication?
figParams = {'units','inches','width',4,'height',2.5,'font','Helvetica','fontsize',10};

%% Plot one: Ptas comparison
ptasFn = 'FINAL_RESULTS/PtasCompare_4.mat';
ptasS  = load(ptasFn, 'problems');
nRange = 4:10;
tRange = 1:12;

bounds(length(nRange), length(tRange)) = 0;

for i = 1:length(nRange)
    subProblems = ptasS.problems{nRange(i)};
    
    [A, B, ~]  = cellfun(@(r) BBP(r.theta, r.W, 0, 0), subProblems, 'UniformOutput', false);
    [~, ~, bound, a, b] = cellfun(@(s, a, b) getIntervalSz(a, b, s.W, 0.01), subProblems, A, B);
    %[~, ~, bound, a, b] = cellfun(@(s, a, b) getIntervalSz(s.mkMisc.A, s.mkMisc.B, s.W, 0.01), subProblems, A, B);

    a > b
    
    bounds(i,:) = bound';              

    rawTimes(i,:) = cellfun(@(p) p.rawMisc.mexTotTime, subProblems);
    bbpTimes(i,:) = cellfun(@(p) p.bbpMisc.mexTotTime, subProblems);
    mkTimes(i,:)  = cellfun(@(p) p.mkMisc.mexTotTime, subProblems);
end

Ns = vec(repmat(nRange, 1, length(tRange)));

% TODO: Use the fig command!
ptasH = fig(figParams{:});
axis([-Inf Inf 5e-1 1e4]);

hold on;
scatter(Ns, 4e-9 * vec(bounds), '+');
scatter(Ns, vec(rawTimes), 'o');
scatter(Ns, vec(bbpTimes), 'x');
scatter(Ns, 20*vec(mkTimes), '^');

title('Runtime vs n');
xlabel('n');
ylabel('Time (s)');
set(gca, 'YScale', 'Log');
lh = legend('Bound', 'PTAS', 'BBP + PTAS', 'MK + PTAS', 'Location', 'SouthEast');
set(lh, 'Box', 'Off');
set(lh, 'FontSize', 8);


export_fig('FINAL_RESULTS/plots/ptas.pdf', ptasH);

%% Plot two: The varying degree, fixed dW plots
nNodes = [8 16 32];
nN     = length(nNodes);

load('FINAL_RESULTS/VaryingRegular_8_RESULTS.mat');
rr8  = [results{:}];
[d8, t8, a8, b8] = varyingTheoryBound(rr8);

load('FINAL_RESULTS/VaryingRegular_16_RESULTS.mat');
rr16 = [results{:}];
[d16, t16, a16, b16] = varyingTheoryBound(rr16);

load('FINAL_RESULTS/VaryingRegular_32_RESULTS.mat');
rr32 = [results{:}];
[d32, t32, a32, b32] = varyingTheoryBound(rr32);

% We need to generate theoretical bounds. We forgot to save A and B, so we
% have to regenerate them.
%

% These files used 1e-4 and 10000 threshold for MK.
%[min8, max8]   = range([rr8.d]);

%[min16, max16] = range([rr16.d]);

maybeH = fig(figParams{:});

hold on;
scatter([rr8.d],  [rr8.mkTime], 16, 'r');
scatter([rr16.d], [rr16.mkTime], 16, 'b');
scatter([rr32.d], [rr32.mkTime], 16, 'g');

legend('N = 8' ,'N = 16', 'N = 32', 'Location', 'Best');

plot(d8, t8, '-.r');
plot(d16, t16, '--.b')
plot(d32, t32, ':.g');

title('Runtime vs degree and N');
xlabel('Degree');
ylabel('Runtime (sec.)');

set(gca, 'YScale', 'Log');

%legend('N = 8' ,'N = 16', 'N =32', 'Location', 'Best');

export_fig('FINAL_RESULTS/plots/maybe.pdf', maybeH);

%% Plot three: Erdos-Renyi random runtimes
%% Setup part:

% Choose a dw 
dw = 8;

nNodeVec = [4 8 16 32 64 128];
%nNodeVec = [8 16 32 64];
textVec  = arrayfun(@(n) ['N = ' num2str(n)], nNodeVec, 'UniformOutput', false);
commaStr = arrayfun(@(n) [num2str(n) ', '], nNodeVec, 'UniformOutput' ,false);
commaStr = [commaStr{:}];
nFiles   = length(nNodeVec);

commaStr(end-1:end) = [];

% Assume uniform across files
nTrials = 100;

kappas(nTrials,nFiles) = 0;
theoryBounds(nTrials,nFiles) = 0;

allMooijTimes(nTrials,nFiles) = 0;
allMaxFlowTimes(nTrials,nFiles) = 0;
allSumMooijGap(nTrials,nFiles) = 0;

totDiffs = cell(nFiles, 1);

for n = 1:nFiles
    %%
    nNodes = nNodeVec(n);
    fn     = sprintf('FINAL_RESULTS/ZeroER/mk_renyi_dw_%d_T_0_nNodes_%d.mat', ...
                     dw, nNodes);
    %fn     = sprintf('FINAL_RESULTS/mk_renyi_const_w_6_T_0_nNodes_%d.mat', ...
    %                 nNodes);    
    load(fn);
    
    % Collate the data we want
    
    % For the "Average error of one-marginals"
    %totDiffs{n} = totDiff;
    
    % Prepare our actual runtimes vs the theoretical bounds
    %%
    for t = 1:nTrials
        % Select the Mooij numbers
        % Remember that epsilon = 0.01
        %[~, kappas(t,n), theoryBounds(t,n)] = getIntervalSz(A(:,t,2),B(:,t,2), problems{t}.W, 0.01);

        % Calculate instead from one BBP iteration
        [A, B, ~] = BBP(problems{t}.theta, problems{t}.W, 0, 0);
        [~, kappas(t,n), theoryBounds(t,n)] = getIntervalSz(A, B, problems{t}.W, 0.01);

        [Amk, Bmk] = BBP_MK_opt_mex(problems{t}.theta, problems{t}.W, 0.01, 1);
        [isz, ~, ~] = getIntervalSz(Amk, Bmk, problems{t}.W, 0.01);
        nIntervals(t,n) = sum(1 - Bmk - Amk) ./ isz;
    end
        
    allMooijTimes(:,n) = betheTimes(:,2);
    allMaxFlowTimes(:,n) = maxFlowTimes(:,2);
    allSumMooijGap(:,n) = sum(ABgap(:,:,2),1)';
end

%save collated.mat allBetheTimes allMooijTimes allSumBetheGap allSumMooijGap kappas theoryBounds;

% %% Plot errors -- IN PUBLICATION AS OF 13 MARCH
% epsilonTxt = [' for \epsilon = ' num2str(epsilon)];
% 
% fig('units','inches','width',6.5,'height',1.75,'font','Helvetica','fontsize',10);
% for n = 1:nFiles
%     subplot(1, nFiles, n);    
%     hist(totDiffs{n}(:,1));
%     title(['Avg. err, n = ' num2str(nNodeVec(n))]);
% end

%% Plot runtime vs nodes -- IN PUBLICATON AS OF 13 MARCH
% fig('units','inches','width',3,'height',3,'font','Helvetica','fontsize',10);
% hold on;
% % TODO: Remove hardcoded labels
% boxplot(allBetheTimes, 'notch', 'marker', 'labels', {'4', '8', '16'}, 'colors', 'g');
% boxplot(theoryBounds, 'notch', 'marker',  'labels', {'', '', ''}, 'colors', 'b');
% set(gca, 'YScale', 'Log');
% title(['Runtime and bound (top) vs nodes; ' num2str(nTrials) ' trials']);
% hold off;

% fig('units','inches','width',3,'height',3,'font','Helvetica','fontsize',12);
% hold on;
% boxplot(allMooijTimes, 'notch', 'marker', 'labels', nodeLabels, 'colors', 'g');
% boxplot(theoryBounds, 'notch', 'marker',  'labels', nodeLabels, 'colors', 'b');
% set(gca, 'YScale', 'Log');
% title(['Mooij runtime (bottom) and theoretical bound (top) vs nodes; ' num2str(nTrials) ' trials']);
% hold off


%% Plot runtime vs Kappa -- IN PUBLICATION AS OF 13 MARCH
boundH = fig(figParams{:});

hold on;
%lineSpecs = {'-.r', '-om', '-xb', '-*g', '-^y', '-vc'};
colors    = {'r', 'c', 'm', 'b', 'g', 'k'};
markers   = {'v', 's', 'o', '*', '^', 'x'};

for n = 1:nFiles        
    [xs, ix] = sort(kappas(:,n));
    ys = allMooijTimes(:,n);
    ys = ys(ix);
    ts = theoryBounds(:,n);
    ts = ts(ix);

    scatter(xs, ys, 16, colors{n}, markers{n});
    set(gca, 'XScale', 'Log');
    set(gca, 'YScale', 'Log');    
 
    lineSpec = ['-' markers{n} colors{n}];
    lh = plot(xs, 1e-4 * ts, lineSpec);
    set(lh, 'MarkerSize', 4);
    xlabel('\Sigma^{3/4}\Omega^{3/2}');
    ylabel('Time (s)');    
end

title(['Runtime vs \Sigma^{3/4}\Omega^{3/2}; n = ' commaStr]);
%title('Runtime vs \kappa');
%legend(textVec{:})
axis([-Inf 1e7 1e-2 Inf]);
set(gca, 'XTick', [1e2 1e4 1e6]);
set(gca, 'YTick', [1e0 1e5 1e10 1e14]);
hold off;

export_fig('FINAL_RESULTS/plots/bound.pdf', boundH);

%% Do the analysis for the random plots
%vmaxFlowTimes = m(:);
vnIntervals = nIntervals(:);
vmaxFlowTimes = allMaxFlowTimes(:);

nzIdxs = nIntervals(:) ~= 0 & vmaxFlowTimes(:) ~= 0;

nodeStats = regstats(log(vmaxFlowTimes(nzIdxs)), log(vnIntervals(nzIdxs)), 'linear');
%edgeStats = regstats(log(vmaxFlowTimes(nzIdxs)), 
