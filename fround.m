function [ y ] = fround( x, d )
% y = fround(x, d) round x to d significant figures.

    y = round(10^d .* x) ./ (10^d);

end

