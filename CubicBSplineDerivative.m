function [out] = CubicBSplineDerivative(input)
%CUBICBSPLINEDERIVATIVE vectorized cubic Bspline slope
%   input dimensions are Nx1


out = NaN(size(input));

c1_ind = find(0 <= abs(input) & abs(input) < 1);   % case1
x = input(c1_ind);
out(c1_ind) = 3/2*x.*abs(x)-2*x;

c2_ind = find(1 <= abs(input) & abs(input) < 2);   % case2
x = input(c2_ind);
out(c2_ind) = -x.*abs(x)/2 + 2*x - 2*x./abs(x);

c3_ind = find(2 <= abs(input));   % case3
out(c3_ind) = 0;
    
% if (0 <= abs(x)) && (abs(x) < 1)
%    	out = 3/2 * x * abs(x) - 2*x;
% elseif (1 <= abs(x)) && (abs(x) < 2)
%     out = -x*abs(x)/2 + 2*x - 2*x/abs(x);
% else
%     out = 0;
% end
end

