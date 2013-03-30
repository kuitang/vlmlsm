function [ J, eta ] = mooijParam(W, theta)
% mooijParam - Convert our paper's W, theta parmeterization to MK's.
%
% [ J, eta ] = mooijParam(W, theta)
%
% NOTE: eta here is theta in MK's paper

% A bit more cumbersome than just inverting Mooij's equations; just write
% out the grid notation
    J = 0.25 * W;
    eta = 0.5*theta + sum(J, 2);
    
end

