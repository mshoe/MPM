function [W_grad] = GridWeightsGradient(Xp, Np, gridDimX, gridDimY, gridX0, gridY0, h)
%GRIDWEIGHTSGRADIENT Summary of this function goes here
%   Detailed explanation goes here

W_grad = NaN(Np, gridDimX, gridDimY, 2);
for p = 1:Np
    W_grad(p,:,:,:) = GridWeightsParticleGradient(Xp(p,:), gridDimX, gridDimY, gridX0, gridY0, h);
end

end

