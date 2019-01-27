function [VOL0p] = InitialParticleVolumes(Wpg, Mp, Np, gridDimX, gridDimY, h)
%INITIALVOLUMES Summary of this function goes here
%   Detailed explanation goes here
Mg = reshape(sum(Wpg.*Mp, 1), [gridDimX gridDimY]);

VOL0p = zeros(Np, 1);
for p = 1:Np
    massFromGrid = sum(sum(squeeze(Wpg(p,:,:)).* Mg));
    if massFromGrid ~= 0.0
        VOL0p(p) = Mp(p) * h^2 / massFromGrid;
    else
        VOL0p(p) = 0.0;
    end
end

end

