function [ y ] = fround( x, d )
% fround Round x to d figures after the decimal point.
%   y = fround(x, d) rounds after the decimal point, NOT literally sigfigs.

    y = round(10^d .* x) ./ (10^d);

end

