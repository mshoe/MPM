function [VOL0p] = InitialParticleVolumes(W, Mp, Np, gridDimX, gridDimY, h)
%INITIALVOLUMES Summary of this function goes here
%   Detailed explanation goes here
Mg = zeros(gridDimX, gridDimY);
for p=1:Np
    Mg = Mg + reshape(W(p,:,:), [gridDimX, gridDimY]) * Mp(p);
end

VOL0p = zeros(Np, 1);
for p = 1:Np
    massFromGrid = sum(sum(squeeze(W(p,:,:)).* Mg));
    if massFromGrid ~= 0.0
        VOL0p(p) = Mp(p) * h^2 / massFromGrid;
    else
        VOL0p(p) = 0.0;
    end
end

end

