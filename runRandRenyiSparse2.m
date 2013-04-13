% JUST RUN IN SERIAL TO AVOID CONTENTION

nNodes = [8 16 32 64];
%nNodes = [64 32 16 8];
N = length(nNodes);

matlabpool open 12

parfor i = 1:N
   % TODO: Run dw = 8 later!

   randRenyiSparse2(nNodes(i), 8, 0);
end

