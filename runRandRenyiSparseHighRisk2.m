% JUST RUN IN SERIAL TO AVOID MEMORY CONTENTION

% BIG IS FASTER THAN DEEP.
%% Second, 128 and 256 node graphs with "reasonable" dw
%nNodesBig = [128 256]
%nBig = length(nNodesBig);
%
%%for i = 1:nBig
%%    for dw = [1 2 4]
%%       fprintf(1, 'ER Graphs; nNodes = %g; dw = %g\n', nNodesBig(i), dw);
%%       randRenyiSparse(nNodesBig(i), dw, 0)
%%   end
%%end
%
%%% Finally, big nodes AND large dw
%for i = 1:nBig
%    for dw = [6 8 10]
%       fprintf(1, 'ER Graphs; nNodes = %g; dw = %g\n', nNodesBig(i), dw);
%       randRenyiSparse(nNodesBig(i), dw, 0);
%   end
%end

%% Not really high-risk, but we forgot
%for dw = [1 2 4]
%    randRenyiSparse(4, dw, 0);
%end
%
%% First: 4 - 16 nodes with dw = 8 and dw = 16
nNodesSmall = [4 8 16 32 64];
nSmall = length(nNodesSmall);

% dw is the bigger limiting factor

for dw = [6 8 10]
    for i = 1:nSmall
       fprintf(1, 'ER Graphs; nNodes = %g; dw = %g\n', nNodesSmall(i), dw);
       randRenyiSparse(nNodesSmall(i), dw, 0);
   end
end

