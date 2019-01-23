function [out] = CubicBSpline(x)
%CUBICBSPLINE Summary of this function goes here
%   Detailed explanation goes here
% inputs:
   % x: value we want to interpolate
x = abs(x);
if (0 <= x) && (x < 1)
    out = 0.5 * x^3 - x^2 + 2/3;
elseif (1 <= x) && (x < 2)
    out = 1/6 * (2 - x)^3;
else
    out = 0;
end
   
end

