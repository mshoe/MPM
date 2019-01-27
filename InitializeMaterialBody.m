function [ Np, Xp, Mp, Vp, Fp, Bp, VOL0p, P_mew, P_lam ] = ...
    InitializeMaterialBody(shape, pMass, pSpacing, youngMod, poisson, V0, gridDimX, gridDimY, gridX0, gridY0, h)
% Each material body needs:
    % Np = # of material points
    % Xp = Positions
    % Mp = Masses
    % Vp = Velocity
    % Fp = Deformation gradient
    % Bp = Affine matrix for APIC transfer
    % VOL0p = Initial volumes
    % P_mew, P_lam = Piola-Kirchoff variables related to Young's Modulus
        % and Poissons Ratio of material
%% 

% generate the point cloud
[Xp, Np] = PointCloud(shape, pSpacing, gridDimX, gridDimY, gridX0, gridY0);    

% mass of material points
Mp = pMass * ones(Np, 1);

% velocity of material points
Vp = repmat(V0, Np, 1);

% deformation gradient at each material point
Fp = NaN(Np, 2, 2);
% originally no deformation:
Fp(:,:,:) = repmat(reshape(eye(2),[1,2,2]), Np, 1, 1);

% initial APIC B matrices
Bp = zeros(Np, 2, 2);

% intial weights (for calculating initial volumes)
Wpg = GridWeights(Xp, Np, gridDimX, gridDimY, h); 

% initial volumes
VOL0p = InitialParticleVolumes(Wpg, Mp, Np, gridDimX, gridDimY, h);

P_mew = youngMod / (2 * (1 + poisson));
P_lam = youngMod * poisson / ((1 + poisson)*(1 - 2*poisson));

end

