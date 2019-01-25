function [W] = GridWeights(Xp, Np, gridDimX, gridDimY, h)
%GRIDWEIGHTS Summary of this function goes here
%   Detailed explanation goes here


W = NaN(Np, gridDimX, gridDimY);
for p = 1:Np
    W(p,:,:) = GridWeightsParticle(Xp(p,:), gridDimX, gridDimY, h);
end
end

