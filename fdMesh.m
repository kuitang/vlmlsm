function [ gams ] = fdMesh(theta, W, A, B, L, U, epsilon, method)
% afdm  First derivative mesh
%   revised by AW, think (hope) I've fixed a few things
%
%   gams = fdMesh(A, B, L, U, epsilon, method) calculates the first-derivative mesh
%   to accuracy epsilon. gams is a N-vector of step sizes, one in each
%   variable. A, B, L, U are N-vectors of bounds; L and U come from BBP.
%   Define \gamma as the *width* of the mesh (so that, e.g. 
%   q_{i,n+1} = q_{i,n} + \gamma_i). (This is at variance with notes2,
%   where \gamma is the distance to the nearest mesh point. So this
%   function will sometimes use 2*\gamma).
%
%   For the moment this means first mesh point is at Ai+gamma_i/2 then keep
%   adding gamma_i to get the next point, until we reach >= Bi
%
%   Note we could easily compute the L and U here, but choose not to for
%   now.
%
%   method can be simple, lagrangian, adaptive. If adaptive, gams will be a
%   length-N cell array, and will contain an array of grid points instead
%   of just step sizes.
%
%   Calculate the maximum first derivative in each direction. The bound
%   is 
%
%   -\theta_i - W_i + \log \frac{q_i}{1 - q_i}
%   \leq \frac{\partial F}{\partial q_i}
%   \leq -\theta_i + V_i + \log \frac{q_i}{1 - q_i}
%
%   To bound, we consider the maximum modulus of the bound on either side.
%   This amounts to checking the lower-left and upper-right corners.
%
%   Later we extend     - integrate rather than sum max
%                       - factor in direction, so e.g. if start from left
%                       Ai, the first mesh point should cover all possible
%                       global min q^ to the left, so only need use
%                       integral of upper bound of derivative (and can
%                       ignore the lower bound function). Then we compute
%                       the 'reach to the right' of the first mesh point,
%                       where for this we must use the lower bound
%                       derivative, which is used for the upper bound of
%                       the - derivative.  From there, again we use the
%                       upper bound only to compute the next mesh point,
%                       which must reach back to the left to the previous
%                       reach point, and continue...

    nNodes = length(theta);

    % We can raise the lower bound curve by \log U_i and lower the upper
    % bound curve by \log \log L_i.    
    Wpos =  sum(W .* (W > 0), 2); % W_i in notes
    Wneg = -sum(W .* (W < 0), 2); % V_i in notes
    
    % Constants; exclude the \log \frac{q_i}{1 - q_i} for now.
    lb   = -theta - Wpos + log(U);
    ub   = -theta + Wneg - log(L);
    
    ll   = lb + log(A) - log(1 - A);
    ur   = ub + log(1 - B) - log(B);
    D    = max(abs(ll), abs(ur));    
    
    S      = 1 - B - A; % gap sizes
    
    % Now, calculate the gamma mesh.
    if strcmp(method, 'simple')
        iszs   = 2 * epsilon / (nNodes * D);  % epsilon and n are scalars, D is a vector
    elseif strcmp(method, 'lagrangian')        
        sqrtSoverD = sqrt(S ./ D); sqrtSD = sqrt(S .* D); sumsqrtSD=sum(sqrtSD);
        iszs   = 2 * epsilon * sqrtSoverD / sumsqrtSD;        
    elseif strcmp(method, 'adaptive')
        % FILL THIS IN
    end
    
    % If we used a non-adaptive method, fill out the mesh    
    if isempty(strfind(method, 'adaptive'))
        gams = cell(nNodes, 1);
        for n = 1:nNodes
            q = A(n):iszs(n):(1 - B(n));
            if isempty(q) % Whole interval is within intervalSz
                q = [A(n) (1 - B(n))];
            elseif q(end) ~= (1 - B(n)) % Did not include the endpoint
                q = [q (1 - B(n))];
            end

            assert(length(q) >= 2, 'not enough intervals!');                
            gams{n} = q;
        end
    end    
    
end
