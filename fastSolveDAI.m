function [ logZ, oneMarginals, logZHist, varargout ] = fastSolveDAI(nNodes, nEdges, psi, method, daiOpts)
% fastSolveDAI Wrap libDAI with marginalization. Must be in psi form.
%   [ logZ, oneMarginals, [twoMarginals] ] = fastSolveDAI(nNodes, nEdges, psi, method, daiOpts)
%
%   psi          - Potentials. Get from makePsi.
%
%   logZ         - exact (JTREE only) or approximate log partition function
%   oneMarginals - nNodes x 1 vector of P(x_n = 1)
%   twoMarginals - 2 x 2 x nEdges array of 2x2 matrices where
%                  M(qi,qj) = P(x_i = qi, x_j = qj).
%                  We only compute pairwise marginals for edges that appear
%                  in W.
 
    % Call DAI
    if strcmp(method, 'BP')
        [logZ,~,~,oneMStruct,twoMStruct,~,logZHist] = dai(psi, method, daiOpts);
    else
        [logZ,~,~,oneMStruct,twoMStruct,~] = dai(psi, method, daiOpts);
        logZHist = 0;
    end        
    
    % oneMStruct is ordered by the nodes because that's how we entered them.
    % At least, we hope so.
    oneMarginals(nNodes,1) = 0;
    for n = 1:nNodes
        %assert(oneMStruct{n}.Member == n);
        oneMarginals(n) = oneMStruct{n}.P(2);
    end
        
    if nargout == 4
        varargout{1}(2,2,nEdges) = 0;
        for ne = 1:nEdges        
            %assert(all(twoMStruct{mStructIdx}.Member == vars(ne,:)));
            varargout{1}(:,:,ne) = twoMStruct{ne + nNodes}.P;
        end
    end
    
end

