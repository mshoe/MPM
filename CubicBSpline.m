function [out] = CubicBSpline(input)
%CUBICBSPLINE vectorized cubic Bspline
%   input dimensions are Nx1

out = NaN(size(input));
input = abs(input);

c1_ind = find(0 <= input & input < 1);   % case1
x = input(c1_ind);
out(c1_ind) = 0.5 * x.^3 - x.^2 + 2/3;

c2_ind = find(1 <= input & input < 2);   % case2
x = input(c2_ind);
out(c2_ind) = 1/6 * (2 - x).^3;

c3_ind = find(2 <= input);   % case3
out(c3_ind) = 0;
   
% x = abs(input);
% if (0 <= x) && (x < 1)
%    	out =  0.5 * x^3 - x^2 + 2/3;
% elseif (1 <= x) && (x < 2)
%     out = 1/6 * (2 - x)^3;
% else
%     out = 0;
% end
end

