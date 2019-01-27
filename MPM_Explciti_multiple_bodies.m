clear all;
clc;

timeSteps = 30;
dt = 1/60;

% number of material points
% Np = 300;



%% Grid initialization

% discretized dimensions
gridDimX = 20;
gridDimY = 20;
h = 1; % grid spacing

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

% mass of material points
pMass = 0.01;
pSpacing = h/4;

% Young's Modulus and Poisson's ratio
youngMod = 1e2;
poisson = 0.20; % cannot be 0.5
V0 = [5.0 0.0];

disp("Generating material body 1");
tic;
c1.x = 12;
c1.y = 8;
c2.x = 15;
c2.y = 10;
r = 2;
shape1 = struct('type', "twoCirclesUnion", 'c1', c1, 'c2', c2, 'r', r);
[ Np, Xp, Mp, Vp, Fp, Bp, VOL0p, P_mew, P_lam ] = ...
    InitializeMaterialBody(shape1, pMass, pSpacing, youngMod, poisson, V0, ...
                            gridDimX, gridDimY, gridX0, gridY0, h);
% Material Point Clouds
N_MPC = 1;
MPC(1) = struct('N', Np, 'X', Xp, 'M', Mp, 'V', Vp, 'F', Fp, 'B', Bp, 'VOL0', VOL0p, 'mew', P_mew, 'lam', P_lam);
toc;

%% Other MPM constants
STICKY = 0.9;
NONSTICKY = 1 - STICKY;

% this matrix doesn't change? APIC D
Dp_i = inv(1/3 * h^2 * eye(2));

%% MPM Algorithm
% ************************************************************** %


pointSize = 4/pSpacing;
fig = figure;%('visible','off');

plot([3 18 18 3 3 18], [3 3 18 18 3 3], 'w', 'LineWidth',7);
set(gca, 'Color', 'k');
ax = gca;
ax.GridColor = [1.0, 1.0, 1.0];
hold on;
props = {'LineStyle','none','Marker','o','MarkerEdge','b','MarkerSize',6};
axis([1 20 1 20]);
xticks(1:20);
yticks(1:20);
grid on;
set(fig, 'Position', [50, 50, 700, 700]);
M = struct('cdata', cell(1,timeSteps+1), 'colormap', cell(1,timeSteps+1));

% Draw the initial configurations of our material point clouds
scatterPoints = gobjects(N_MPC, 1);
for m = 1:N_MPC
    scatterPoints(m) = scatter(Xp(:, 1), Xp(:,2), pointSize, 'filled', 'red');
end
M(1) = getframe;



for t=1:timeSteps
for m = 1:N_MPC
        t
    frameStart = tic;
    %% 0. Compute grid weights
    % first get the grid weights for each material point
    Wpg = GridWeights(Xp, Np, gridDimX, gridDimY, h);
    Wpg_grad = GridWeightsGradient(Xp, Np, gridDimX, gridDimY, h);

    %% 1. Particle to grid transfer (using APIC)
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

    %% 2. Compute grid velocities
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
    %% 3. Identify grid degrees of freedom
    [degX, degY] = find(Mg);
    GridDegrees = [degX, degY];
    numDegrees = size(degX);
    numDegrees = numDegrees(1);
    GInds1 = sub2ind([gridDimX gridDimY 2], degX, degY, ones(numDegrees,1));
    GInds2 = sub2ind([gridDimX gridDimY 2], degX, degY, 2*ones(numDegrees,1));
    %% 4. Compute explicit grid forces
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
    end
    %% 5. Grid velocity update

    Vg(GInds1) = Vg(GInds1) + dt*force_g(GInds1) ./ Mg(GInds1);
    Vg(GInds2) = Vg(GInds2) + dt*force_g(GInds2) ./ Mg(GInds1) - dt*9.8;
end
for m=1:N_MPC
    %% 5.5 Compute collisions
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
    %% 6. Update particle deformation gradient
    for p=1:Np
        PGInds1 = sub2ind([Np gridDimX gridDimY 2], p*ones(numDegrees,1), degX, degY, 1*ones(numDegrees,1));
        PGInds2 = sub2ind([Np gridDimX gridDimY 2], p*ones(numDegrees,1), degX, degY, 2*ones(numDegrees,1));

        temp = NaN(2,2);
        temp(1,1) = sum(Vg(GInds1) .* Wpg_grad(PGInds1));
        temp(1,2) = sum(Vg(GInds1) .* Wpg_grad(PGInds2));
        temp(2,1) = sum(Vg(GInds2) .* Wpg_grad(PGInds1));
        temp(2,2) = sum(Vg(GInds2) .* Wpg_grad(PGInds2));

          Fp(p,:,:) = (eye(2) + dt * temp) * reshape(Fp(p,:,:), [2 2 1]);
    end

    %% 7. Grid to particle transfer
    for p=1:Np

        PGInds = sub2ind([Np gridDimX gridDimY], p*ones(numDegrees,1), degX, degY);

        Vp(p,1) = sum(Wpg(PGInds).*Vg(GInds1));
        Vp(p,2) = sum(Wpg(PGInds).*Vg(GInds2));

        dXgp1 = Xg(GInds1) - Xp(p,1);
        dXgp2 = Xg(GInds2) - Xp(p,2);
        Bp(p,1,1) = sum(Wpg(PGInds).*Vg(GInds1).*dXgp1);
        Bp(p,1,2) = sum(Wpg(PGInds).*Vg(GInds1).*dXgp2);
        Bp(p,2,1) = sum(Wpg(PGInds).*Vg(GInds2).*dXgp1);
        Bp(p,2,2) = sum(Wpg(PGInds).*Vg(GInds2).*dXgp2);
    end
    %% 8. Particle advection
    Xp = Xp + dt*Vp;
end
frameTime = toc(frameStart);
fprintf("MPM timestep took %f seconds.\n", frameTime);

%% Plot/Animation
tic;

% draw the updated material point clouds
for m = N_MPC
    delete(scatterPoints(m));
    scatterPoints(m) = scatter(Xp(:, 1), Xp(:,2), pointSize, 'filled', 'red');
end

M(t+1) = getframe;
drawTime = toc;
fprintf("Drawing point clouds took %f seconds.\n", drawTime);
end

fprintf("\n\nMPM complete, now writing video...\n");
tic;
v = VideoWriter('mpm_test.avi');
v.FrameRate = 60;
open(v);
writeVideo(v, M);
close(v);
toc;