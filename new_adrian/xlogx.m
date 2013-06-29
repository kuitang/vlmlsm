function [ f ] = xlogx( x )
%[ f ] = xlogx( x )
%   Returns x * log x, using natural log and returns 0 if x=0
if x==0
    f=0;
else
    f=x*log(x);
end
end

