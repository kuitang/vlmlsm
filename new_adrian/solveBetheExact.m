function [ mu, xi, energy ] = solveBetheExact(theta, W, varargin)
% [mu, xi, energy] = solveBetheExact(theta, W, epsilon)
%
%   Find (global) minimum of the Bethe free energy.
%
%   Right now, parameterizes an ILP using our bounds. But this will seem
%   silly. We should use gradient information. What about alpha-BB?
%
%   NOTE: Energy is -log Z. But the outer (dual decomposition) loop
%   introduces some constant shift, which it takes care of.

    % FILL IN THE DETAILS; 

    p = inputParser;
    p.addRequired('theta', @isnumeric);
    p.addRequired('W');    
    p.addParamValue('epsilon', 0.1); % It's gonna be slow...
    p.addParamValue('algo', 'gurobi');
    p.addParamValue('method', 'adaptiveminsum');    
    p.parse(theta, W, varargin{:});
    
    maxIter = 1000;
    
    o = p.Results;
    N = length(theta);
    
    W = triu(W) + triu(W)';    
    
    if strcmp(o.algo, 'gurobi')
        % Directly solve the ILP with Gurobi
        
        % Get the gams from somewhere...
        [A, B] = bpbound(N, theta, W, maxIter);        
        [gams, nPts] = fdm(theta, W, A, B, o.epsilon, o.method);
        
        % Check the bounds        
        
        [D, ~, Vi, Vm] = boundMRFNew(theta, W, gams);
        L = cellfun(@length, D);        
        
        % NOTE: The newW argument is not used. We need to clean up the code
        % to not really rely on the MultiLabelSubmodular crap... it's kinda
        % dirty right now.
        
        % We add the trivial (variable and edge) constraints; see
        % [Komodakis11], expressed in the Ax = b. Note that x is a BINARY
        % vector; the label for a given node is expressed as a 1-hot
        % subvector of x.
        
        % In particular, each label for each variable has its index in x,
        % as does each pairwise label. Instead of trying to calculate
        % indices, we'll just precompute them.
        %
        % iTheta{u}(k) is the index of the the kth label on node n.
        % iW{u,v}(k,l) is the index of (k,l) label on edge (u,v).
        
        iTheta = cell(N, 1);
        obj    = [];
        
        start = 1;
        for u = 1:N            
            iTheta{u} = start:(start + L(u) - 1);
            start = iTheta{u}(end) + 1;
            
            % Set the objective coefficients (which are just the singleton
            % potentials)
            obj(iTheta{u}) = D{u};
        end
        
        % iW is UT.
        iW = cell(N, N);
        for u = 1:N
            for v = (u+1):N                                
                nLabels = L(u) * L(v);                
                seq     = start:(start + nLabels - 1);
                iW{u,v} = reshape(seq, L(u), L(v));                                
                
                start = seq(end) + 1;
                
                % Set the objective coefficients (which are just the
                % pairwise potentials)
                %
                % TODO: Extract this; put it in a better place... move all
                % this LP stuff to a different function so its easier to
                % debug (you can enter some Toy problems)
                obj(seq) = Vm{Vi(u,v)}(:);
            end
        end
        
        % Objective (each entry is just the value of its potential)        
        
        % WARNING: THIS IS INEFFICIENT. (But the integer program will be
        % way more inefficient.)
        A = sparse(0);
        rhs = [];
        
        iConstr = 1;        
        % Single label constraints
        for u = 1:N                        
            for k = 1:L(u)
                % The constraint (sum of the subvector of x for node n
                % equals 1)                
                A(iConstr, iTheta{u}) = 1;                
                rhs(iConstr) = 1;
                
                iConstr = iConstr + 1;
            end
        end
        
        for u = 1:N
            for v = (u + 1):N
                % Row-wise marginalization constraints                
                for ku = 1:L(u)
                    % Fix label ku and marginalize v
                    A(iConstr, iW{u,v}(ku,:)) = 1;
                    
                    % Set equality to ku label of u
                    A(iConstr, iTheta{u}(ku)) = -1;
                    rhs(iConstr) = 0;
                    
                    iConstr = iConstr + 1;                    
                end
                
                % Column-wise marginalization constraints
                for kv = 1:L(v)
                    % Fix label kv and marginalize u
                    A(iConstr, iW{u,v}(:,kv)) = 1;
                    
                    % Set equality to the kv label of v
                    A(iConstr, iTheta{v}(kv)) = -1;
                    rhs(iConstr) = 0;
                    
                    iConstr = iConstr + 1;
                end
            end
        end
        
        % CHECK EVERYTHING HERE
        %A
        %rhs
        %obj
        
        % Now do the optimization
        sense = repmat('=', 1, iConstr - 1);
        vtype = repmat('B', 1, length(obj));
        
        model = var2struct(A, obj, rhs, vtype, sense);
        params = struct('OutputFlag', 0);
        
        result = gurobi(model, params);
        assert(strcmp(result.status, 'OPTIMAL'), 'ILP not solved to optimality!');
        energy = result.objval;
        
        % Recover the (integer) labels
        label(N) = 0;
        for u = 1:N
            k = find(result.x(iTheta{u}));
            assert(length(k) == 1, 'Single label constraint failed!');
            label(u) = k;
        end
        
        % Recover energies
        mu(N) = 0;
        xi = zeros(N,N);
        
        
        % XI IS NOT VM!! There is a formula. (In fact, closed form; it's
        % not even a silly energy calculation)
        iEdge = 1;
        alpha = exp(W) - 1;
        for u = 1:N
            assert(length(gams{u}) == L(u) && length(D{u}) == L(u), 'marginal, energy, and label mismatch');
            mu(u) = gams{u}(label(u));
        end
        
        for u = 1:N
            for v = (u+1):N
                [~, xi(u,v)] = marginalize(alpha(u,v), mu(u), mu(v));
            end
        end
    else
        error('solveBetheExact:algo', 'Algo %s not recognized', o.algo);
    end

end

