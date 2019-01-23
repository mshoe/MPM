function [W] = GridWeightsParticle(Xp, gridDimX, gridDimY, gridX0, gridY0, h)
%GRIDWEIGHTS Calculates the weights wij for the particle p
%   Xp = [x y] position of particle

N = @CubicBSpline;

xp = Xp(1);
yp = Xp(2);

W = NaN(gridDimX, gridDimY);
%Xij = repmat((gridX0:1:gridDimX)', 1, gridDimY);
%Yij = repmat((gridY0:1:gridDimY), gridDimX, 1);

for i=1:gridDimX
    for j=1:gridDimY
        xij = gridX0 + i - 1;
        yij = gridY0 + j - 1;
        W(i,j) = N(1/h*(xp - xij)) * N(1/h*(yp - yij));
    end
end
end

