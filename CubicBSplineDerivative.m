function [out] = CubicBSplineDerivative(x)
%CUBICBSPLINEDERIVATIVE Summary of this function goes here
%   Detailed explanation goes here

if (0 <= abs(x)) && (abs(x) < 1)
   	out = 3/2 * x * abs(x) - 2*x;
elseif (1 <= abs(x)) && (abs(x) < 2)
    out = -x*abs(x)/2 + 2*x - 2*x/abs(x);
else
    out = 0;
end
end

