% JUST RUN IN SERIAL TO AVOID CONTENTION

nNodes = [8 16 32 64];
N = length(nNodes);

for i = 1:N
   % TODO: Run dw = 8 later!
   for dw = [1 2 4]
       fprintf(1, 'ER Graphs; nNodes = %g; dw = %g\n', nNodes(i), dw);
       randRenyiSparse(nNodes(i), dw, 0);
   end
end

