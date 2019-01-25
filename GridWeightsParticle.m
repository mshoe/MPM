function [W] = GridWeightsParticle(Xp, gridDimX, gridDimY, h)
%GRIDWEIGHTS Calculates the weights wij for the particle p
%   Xp = [x y] position of particle

N = @CubicBSpline;

xp = Xp(1);
yp = Xp(2);

W = NaN(gridDimX, gridDimY);
[IndX, IndY] = find(ones(gridDimX, gridDimY));
Inds = squeeze(1:1:(gridDimX*gridDimY));
W(Inds) = N(1/h*(xp - IndX)) .* N(1/h*(yp - IndY));
% for i=1:gridDimX
%     for j=1:gridDimY
%         xij = i;
%         yij = j;
%         W(i,j) = N(1/h*(xp - xij)) * N(1/h*(yp - yij));
%     end
% end
end

