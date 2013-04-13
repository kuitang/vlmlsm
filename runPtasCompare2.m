%% Global (fixed) variables 

% Set dw so hard that the problems take nontrivial time.
dw = 6;

% Load the graphs etc.
fn       = 'PtasCompare_4.mat'
problems = cell(10);

matlabpool open 12

for nNodes = 4:10
    b = 2 * floor(log(nNodes));
    infn = sprintf('bregular%d.%d.mat', nNodes, b);
    load(infn);
    nGraphs = length(graphs);

    problems{nNodes} = parPtasCompare(graphs, dw / b);

    save(fn);
end

