function [ mrf ] = fdReduction(theta, W, varargin)
% fdReduction  Construct factor graph to solve binary MRF
%   
%   [ mrf ] = fdReduction(theta, W, param1, value1, ...) blah blah
%
%   mrf is a MIN-SUM problem, a structure with fields D, newW, Vi, Vm, gams
%   (documented below)

    p = inputParser;
    p.addRequired('theta');
    p.addRequired('W');
    p.addParamValue('epsilon', 0.01);
    p.addParamValue('bbpThresh', 1e-4);
    p.addParamValue('bbpMaxIter', 400);
    p.addParamValue('bbpMethod', 'WJ'); % can also be 'MK'
    p.addParamValue('meshMethod', 'simple');
    p.parse(theta, W, varargin{:});
    
    o = p.Results;
    
    nNodes = length(theta);
    
    if strcmp(o.bbpMethod, 'WJ')
        [A, B, alpha, L, U] = BBPNew(theta, W, 'thresh', o.bbpThresh, ...
            'maxIter', o.bbpMaxIter);
    elseif strcmp(o.bbpMethod, 'MK')
        [Ainit, Binit] = bpbound(nNodes, theta, W, o.bbpMaxIter);
        [A, B, alpha, L, U] = BBPNew(theta, W, 'thresh', o.bbpThresh, ...
            'maxIter', o.bbpMaxIter, 'A', Ainit, 'B', Binit);
    else
        error('fdReduction:bbpMethod', 'Unrecognized bbpMethod %s', o.bbpMethod);
    end
    
    gams = fdMesh(theta, W, A, B, L, U, o.epsilon, o.meshMethod);
    
    % D is an N-cell of vectors with D{n}(i) = energy of gams{n}(i)
    %
    % newW is a sparse matrix for the new MRF
    %
    % Vi is sparse matrix of the same dimensions as newW; Vi(i,j) denotes
    % the linear index into Vm where the potential for edge (i,j) is
    % stored.
    %
    % Vm is a cell vector of matrices. Rows (dimension 1) index node i and
    % columns (dimension 2) index node j.
    [D, newW, Vi, Vm, gams] = boundMRFNew(theta, W, gams);
    mrf = var2struct(D, newW, Vi, Vm, gams);
end

