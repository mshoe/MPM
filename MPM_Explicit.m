clear all;
clc;

timeSteps = 600;
dt = 1/60;

% number of material points
% Np = 300;

%% Grid initialization

% discretized dimensions
gridDimX = 20;
gridDimY = 20;
h = 1; % grid spacing
d = h/4;
gridX0 = 1;
gridY0 = 1;
% this means nodes are located at (1,1), (1,2), ... , (20, 20)

BORDER_MIN = 3;
BORDER_MAX = 18;

% position of grid nodes (doesn't change)
Xg = zeros(gridDimX, gridDimY, 2);
for i=1:gridDimX
    for j=1:gridDimY
        Xg(i,j, :) = [gridX0 + i - 1, gridY0 + j - 1];
    end
end

% mass of grid nodes
Mg = zeros(gridDimX, gridDimY);

% momentum of grid nodes
MOMg = zeros(gridDimX, gridDimY, 2);

% Grid velocities
Vg = zeros(gridDimX, gridDimY, 2);

%% Initialization of Material Points


% reference configuration of material points
%X0p = zeros(Np, 2);%, 'gpuArray');
% BoxBounds = [5 8];
% X0p = BoxBounds(1) + rand(Np, 2) * (BoxBounds(2) - BoxBounds(1));
tic;
[X0p, Np] = PointCloud(d, gridDimX, gridDimY, gridX0, gridY0);
toc;

% current config of material points
Xp = X0p;

% mass of material points
pointMass = 0.01;
Mp = pointMass * ones(Np, 1);%, 'gpuArray');

% velocity of material points
Vp = zeros(Np, 2);%, 'gpuArray');

% deformation gradient at each material point
Fp = NaN(Np, 2, 2);%, 'gpuArray');
% originally no deformation:
Fp(:,:,:) = repmat(reshape(eye(2),[1,2,2]), Np, 1, 1);
% originally stretched
% for p=1:Np
%     Fp(p,:,:) = 2*eye(2);
% end

% Young's Modulus and Poisson's ratio
youngMod = 1e2;
poisson = 0.20; % cannot be 0.5
P_mew = youngMod / (2 * (1 + poisson));
P_lam = youngMod * poisson / ((1 + poisson)*(1 - 2*poisson));

STICKY = 0.9;
NONSTICKY = 1 - STICKY;
%% MPM final initializations

% intial weights
Wpg = GridWeights(X0p, Np, gridDimX, gridDimY, h); 

% initial volumes
VOL0p = InitialParticleVolumes(Wpg, Mp, Np, gridDimX, gridDimY, h);

% initial APIC B matrices
Bp = zeros(Np, 2, 2);

% this matrix doesn't change? APIC D
Dp_i = inv(1/3 * h^2 * eye(2));

%% MPM Algorithm
% ************************************************************** %


pointSize = 4/d;
fig = figure;%('visible','off');

plot([3 18 18 3 3 18], [3 3 18 18 3 3], 'w', 'LineWidth',7);
set(gca, 'Color', 'k');
ax = gca;
ax.GridColor = [1.0, 1.0, 1.0];
hold on;
props = {'LineStyle','none','Marker','o','MarkerEdge','b','MarkerSize',6};
MPs = scatter(X0p(:, 1), X0p(:,2), pointSize, 'filled', 'red');
%MPs = line([X0p(:,1), X0p(:,1)], [X0p(:,2), X0p(:,2)], props{:});
axis([1 20 1 20]);
xticks(1:20);
yticks(1:20);
grid on;
set(fig, 'Position', [50, 50, 700, 700]);
M = struct('cdata', cell(1,timeSteps+1), 'colormap', cell(1,timeSteps+1));
M(1) = getframe;



for t=1:timeSteps
    t
frameStart = tic;
%% 0. Compute grid weights
tic;
% first get the grid weights for each material point
Wpg = GridWeights(Xp, Np, gridDimX, gridDimY, h);
Wpg_grad = GridWeightsGradient(Xp, Np, gridDimX, gridDimY, h);
toc;

%% 1. Particle to grid transfer (using APIC)
tic;
% then transfer the masses of the mps to the grid
Mg = reshape(sum(Wpg.*Mp, 1), [gridDimX gridDimY]);

% transfer the momentum of the mps to the grid using APIC formula
MOMg = zeros(gridDimX, gridDimY, 2);
for p=1:Np
    
    wgp_mp = reshape((Wpg(p,:,:)) * Mp(p), [gridDimX gridDimY]);
    dxg_xp = Xg - reshape(Xp(p,:), [1,1,2]);
    BpDp_i = reshape(Bp(p,:,:), [2, 2]) * Dp_i;
    
    BpDp_i = reshape(BpDp_i, [1, 1, 2, 2]);
    
    temp = zeros(gridDimX, gridDimY, 2);
    temp(:,:,1) = BpDp_i(:,:,1,1).*dxg_xp(:,:,1) + BpDp_i(:,:,1,2).*dxg_xp(:,:,2);
    temp(:,:,2) = BpDp_i(:,:,2,1).*dxg_xp(:,:,1) + BpDp_i(:,:,2,2).*dxg_xp(:,:,2);
    MOMg = MOMg + wgp_mp .* (reshape(Vp(p,:), [1,1,2]) + temp);
end
toc;

%% 2. Compute grid velocities


tic;
% non-vectorized, need to find vectorized solution
for i=1:gridDimX
    for j=1:gridDimY
        if Mg(i,j) == 0.0
            Vg(i,j,:) = 0;
        else
            Vg(i,j,:) = MOMg(i,j,:) / Mg(i,j);
        end
    end
end
toc;
%% 3. Identify grid degrees of freedom
tic;
[degX, degY] = find(Mg);
GridDegrees = [degX, degY];
numDegrees = size(degX);
numDegrees = numDegrees(1);
GInds1 = sub2ind([gridDimX gridDimY 2], degX, degY, ones(numDegrees,1));
GInds2 = sub2ind([gridDimX gridDimY 2], degX, degY, 2*ones(numDegrees,1));
toc;
%% 4. Compute explicit grid forces
tic;
% compute determinant
Jp = Fp(:,1,1).*Fp(:,2,2) - Fp(:,1,2).*Fp(:,2,1);

% compute F inverse transpose
Fp_it = NaN(Np, 2, 2);
for p=1:Np
    Fp_it(p,:,:) = inv(reshape(Fp(p,:,:), [2 2]))'; 
end

% compute F transpose
Fp_t = NaN(Np, 2, 2);
for p=1:Np
    Fp_t(p,:,:) = reshape(Fp(p,:,:), [2 2])';
end

% compute Piola Kirchoff stress
Pp = P_mew * (Fp - Fp_it) + P_lam * log(Jp).*Fp_it;

force_g = zeros(gridDimX, gridDimY, 2);
for k=1:numDegrees
    i = GridDegrees(k,1);
    j = GridDegrees(k,2);
    
    % element wise matrix*matrix multiplication
    temp1 = NaN(Np,2,2);
    temp1(:,1,1) = Pp(:,1,1).*Fp_t(:,1,1) + Pp(:,1,2).*Fp_t(:,2,1);
    temp1(:,1,2) = Pp(:,1,1).*Fp_t(:,1,2) + Pp(:,1,2).*Fp_t(:,2,2);
    temp1(:,2,1) = Pp(:,2,1).*Fp_t(:,1,1) + Pp(:,2,2).*Fp_t(:,2,1);
    temp1(:,2,2) = Pp(:,2,1).*Fp_t(:,1,2) + Pp(:,2,2).*Fp_t(:,2,2);
    
    % element wise matrix * vector multiplication
    temp2 = NaN(Np,2);
    temp2(:,1) = temp1(:,1,1).*Wpg_grad(:,i,j,1) + temp1(:,1,2).*Wpg_grad(:,i,j,2);
    temp2(:,2) = temp1(:,2,1).*Wpg_grad(:,i,j,1) + temp1(:,2,2).*Wpg_grad(:,i,j,2);
    
    force_g(i,j,:) = -sum(VOL0p.*temp2,1);
%     for p=1:Np
%         temp = squeeze(Pp(p,:,:)) * squeeze(Fp_t(p,:,:)) * squeeze(Wpg_grad(p,i,j,:));
%         force_g(i,j,:) = force_g(i,j,:) - reshape(VOL0p(p)*temp, [1 1 2]);
%     end
end
toc;
%% 5. Grid velocity update
tic;
% for k=1:numDegrees
%     i = GridDegrees(k,1);
%     j = GridDegrees(k,2);
%     Vg(i,j,:) = Vg(i,j,:) + dt*force_g(i,j,:) / Mg(i,j);
%     Vg(i,j,2) = Vg(i,j,2) - dt * 9.8;
% end

Vg(GInds1) = Vg(GInds1) + dt*force_g(GInds1) ./ Mg(GInds1);
Vg(GInds2) = Vg(GInds2) + dt*force_g(GInds2) ./ Mg(GInds1) - dt*9.8;
%Vg(degX, degY, :) = dt*force_g(degX, degY, :) ./ Mg(degX, degY) - reshape([0, dt*9.8], [1 1 2]);

% compute grid node to border collisions
for k=1:numDegrees
    i = GridDegrees(k,1);
    j = GridDegrees(k,2);
    
    newNodePos = [i*h; j*h] + dt * reshape(Vg(i,j,:), [2 1 1]);
    
    % left right borders
    if (newNodePos(1) < BORDER_MIN) || (newNodePos(1) > BORDER_MAX)
        Vg(i, j, 1) = 0;
        Vg(i, j, 2) = Vg(i, j, 2) * NONSTICKY;
    end
    
    % top bottom borders
    if (newNodePos(2) < BORDER_MIN) || (newNodePos(2) > BORDER_MAX)
        Vg(i, j, 1) = Vg(i, j, 1) * NONSTICKY;
        Vg(i, j, 2) = 0;
    end
end
toc;
%% 6. Update particle deformation gradient
tic;
for p=1:Np
%     temp = zeros(2,2);
%     for k=1:numDegrees
%         i = GridDegrees(k,1);
%         j = GridDegrees(k,2);
% 
%         temp = temp + reshape(Vg(i,j,:), [2 1 1]) * reshape(Wpg_grad(p,i,j,:), [2 1 1 1])';
%     end
    PGInds1 = sub2ind([Np gridDimX gridDimY 2], p*ones(numDegrees,1), degX, degY, 1*ones(numDegrees,1));
    PGInds2 = sub2ind([Np gridDimX gridDimY 2], p*ones(numDegrees,1), degX, degY, 2*ones(numDegrees,1));
    
    temp = NaN(2,2);
    temp(1,1) = sum(Vg(GInds1) .* Wpg_grad(PGInds1));
    temp(1,2) = sum(Vg(GInds1) .* Wpg_grad(PGInds2));
    temp(2,1) = sum(Vg(GInds2) .* Wpg_grad(PGInds1));
    temp(2,2) = sum(Vg(GInds2) .* Wpg_grad(PGInds2));
    
      Fp(p,:,:) = (eye(2) + dt * temp) * reshape(Fp(p,:,:), [2 2 1]);
end
toc;

%% 7. Grid to particle transfer
tic;
for p=1:Np
%     temp = [0 0];
%     for k=1:numDegrees
%         i = GridDegrees(k,1);
%         j = GridDegrees(k,2);
%         
%         temp = temp + Wpg(p,i,j) * reshape(Vg(i,j,:), [2 1 1])';
%     end
%     Vp(p,:) = temp;
    PGInds = sub2ind([Np gridDimX gridDimY], p*ones(numDegrees,1), degX, degY);
    
    Vp(p,1) = sum(Wpg(PGInds).*Vg(GInds1));
    Vp(p,2) = sum(Wpg(PGInds).*Vg(GInds2));

%     temp = zeros(2);
%     for k=1:numDegrees
%         i = GridDegrees(k,1);
%         j = GridDegrees(k,2);
%         
%         temp = temp + Wpg(p,i,j) * reshape(Vg(i,j,:), [2 1 1]) * reshape(Xg(i,j,:) - reshape(Xp(p,:), [1 1 2]), [2 1 1])';
%     end
%     Bp(p,:,:) = temp;
    dXgp1 = Xg(GInds1) - Xp(p,1);
    dXgp2 = Xg(GInds2) - Xp(p,2);
    Bp(p,1,1) = sum(Wpg(PGInds).*Vg(GInds1).*dXgp1);
    Bp(p,1,2) = sum(Wpg(PGInds).*Vg(GInds1).*dXgp2);
    Bp(p,2,1) = sum(Wpg(PGInds).*Vg(GInds2).*dXgp1);
    Bp(p,2,2) = sum(Wpg(PGInds).*Vg(GInds2).*dXgp2);
end
toc;
%% 8. Particle advection
tic;
Xp = Xp + dt*Vp;
toc;
%% Plot/Animation
tic;

delete(MPs);
MPs = scatter(Xp(:, 1), Xp(:,2), pointSize, 'filled', 'red');
%MPs = line([Xp(:,1), Xp(:,1)], [Xp(:,2), Xp(:,2)], props{:});
M(t+1) = getframe;
toc;

toc(frameStart)
end

disp("MPM complete, now writing video...");
tic;
v = VideoWriter('mpm_test.avi');
v.FrameRate = 60;
open(v);
writeVideo(v, M);
close(v);
toc;