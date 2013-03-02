function V = makeMonge(N, M, figs)
% makeMonge Generate random Monge matrix.
%
%   V = makeMonge(N, M, sigFigs) makes a submodular matrix of size NxM with
%   figs digits after the decimal place, if specified.
%
%   
%   We simply generate a random nonnegative "density matrix" and sum it to
%   get a "distribution" Monge matrix.
%
%   References
%   [1] Rainer E. Burkard, "Monge properties, discrete convexity and
%       applications", European Journal of Operational Research, 176: 1
%       (2007), pp 1-14, doi:10.1016/j.ejor.2005.04.050


    d = rand(N, M);
    
    if nargin > 2 && figs
        d = fround(d, figs);
    end
    
    V = zeros(N, M);
    for j = 1:M
        for i = 1:N
            % Burkard (3); note the summation just means "add everything in
            % the lower-left corner"
            V(i,j) = sum(vec(d(i:N,1:j)));
        end
    end
        
    assert(IsMonge(V));
end
