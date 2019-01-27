clear all;
clc;

timeSteps = 600;
dt = 1/60;

% number of material points
% Np = 300;



%% Grid initialization

% discretized dimensions
gridDimX = 40;
gridDimY = 40;
h = 1; % grid spacing

gridX0 = 1;
gridY0 = 1;
% this means nodes are located at (1,1), (1,2), ... , (20, 20)

BORDER_MIN_X = gridX0 + 2;
BORDER_MAX_X = gridDimX - 2;
BORDER_MIN_Y = gridY0 + 2;
BORDER_MAX_Y = gridDimY - 2;

% position of grid nodes (doesn't change)
Xg = zeros(gridDimX, gridDimY, 2);
for i=1:gridDimX
    for j=1:gridDimY
        Xg(i,j, :) = [gridX0 + i*h - 1, gridY0 + j*h - 1];
    end
end

% mass of grid nodes
Mg = zeros(gridDimX, gridDimY);

% momentum of grid nodes
MOMg = zeros(gridDimX, gridDimY, 2);

% Grid velocities
Vg = zeros(gridDimX, gridDimY, 2);

% Grid forces
force_g = zeros(gridDimX, gridDimY, 2);

init_grid = struct('M', Mg, 'V', Vg, 'MOM', MOMg, 'force', force_g, ...
                    'numDegrees', 0, 'degX', [], 'degY', [], 'GridDegrees', [], 'GInds1', [], 'GInds2', []);

%% Initialization of Material Points

% mass of material points
pMass = 0.01;
pSpacing = h/4;

% Young's Modulus and Poisson's ratio
youngMod = 1e2;
poisson = 0.20; % cannot be 0.5



% Material Point Clouds
N_MPC = 2;

disp("Generating material body 1");
tic;
c1.x = 12.5+10;
c1.y = 8+5;
c2.x = 15+10;
c2.y = 10+5;
r = 3;
V0 = [-5.0 2.0];
shape1 = struct('type', "twoCirclesUnion", 'c1', c1, 'c2', c2, 'r', r);
[ Np, Xp, Mp, Vp, Fp, Bp, VOL0p, Wpg, Wpg_grad, P_mew, P_lam ] = ...
    InitializeMaterialBody(shape1, pMass, pSpacing, youngMod, poisson, V0, ...
                            gridDimX, gridDimY, gridX0, gridY0, h);
MPC(2) = struct('N', Np, 'X', Xp, 'M', Mp, 'V', Vp, 'F', Fp, 'B', Bp, 'VOL0', VOL0p, ...
                'Wpg', Wpg, 'Wpg_grad', Wpg_grad, ...
                'mew', P_mew, 'lam', P_lam, ...
                'Color', 'r');
MG(2) = init_grid;
toc;

disp("Generating material body 2");
tic;
c.x = 7;
c.y = 12;
r1 = 1;
r2 = 5;%3;
V0 = [5.0 3.0];
shape2 = struct('type', "demoDonut", 'c', c, 'r1', r1, 'r2', r2);
[ Np, Xp, Mp, Vp, Fp, Bp, VOL0p, Wpg, Wpgr_grad, P_mew, P_lam ] = ...
    InitializeMaterialBody(shape2, pMass, pSpacing, youngMod, poisson, V0, ...
                            gridDimX, gridDimY, gridX0, gridY0, h);
MPC(1) = struct('N', Np, 'X', Xp, 'M', Mp, 'V', Vp, 'F', Fp, 'B', Bp, 'VOL0', VOL0p, ...
                'Wpg', Wpg, 'Wpg_grad', Wpg_grad, ...
                'mew', P_mew, 'lam', P_lam, ...
                'Color', 'b');
MG(1) = init_grid;
toc;

%% Other MPM constants
STICKY = 0.9;
NONSTICKY = 1 - STICKY;

% this matrix doesn't change? APIC D
Dp_i = inv(1/3 * h^2 * eye(2));

%% MPM Algorithm
% ************************************************************** %


pointSize = 4/pSpacing * 40/(gridDimX+gridDimY);
fig = figure;%('visible','off');

%plot([BORDER_MIN_X BORDER_MAX_X BORDER_MAX_X BORDER_MIN_X BORDER_MIN_X BORDER_MAX_X], ...
%    [BORDER_MIN_Y BORDER_MIN_Y BORDER_MAX_Y BORDER_MAX_Y BORDER_MIN_Y BORDER_MIN_Y], 'w', 'LineWidth',7);
set(gca, 'Color', 'k');
ax = gca;
ax.GridColor = [1.0, 1.0, 1.0];
hold on;
axis([gridX0 gridDimX gridY0 gridDimY]);
xticks(gridX0:gridDimX);
yticks(gridY0:gridDimY);
grid on;
set(fig, 'Position', [50, 50, 700, 700]);
M = struct('cdata', cell(1,timeSteps+1), 'colormap', cell(1,timeSteps+1));

% Draw the initial configurations of our material point clouds
scatterPoints = gobjects(N_MPC, 1);
for m = 1:N_MPC
    scatterPoints(m) = scatter(MPC(m).X(:, 1), MPC(m).X(:,2), pointSize, 'filled', MPC(m).Color);
end
M(1) = getframe;



for t=1:timeSteps
fprintf("Timestep t = %d.\n", t);
frameStart = tic;
for m = 1:N_MPC
    %% 0. Compute grid weights
    % first get the grid weights for each material point
    MPC(m).Wpg = GridWeights(MPC(m).X, MPC(m).N, gridDimX, gridDimY, h);
    MPC(m).Wpg_grad = GridWeightsGradient(MPC(m).X, MPC(m).N, gridDimX, gridDimY, h);

    %% 1. Particle to grid transfer (using APIC)
    % then transfer the masses of the mps to the grid
    MG(m).M = reshape(sum(MPC(m).Wpg.*MPC(m).M, 1), [gridDimX gridDimY]);

    % transfer the momentum of the mps to the grid using APIC formula
    MG(m).MOM = zeros(gridDimX, gridDimY, 2);
    for p=1:MPC(m).N

        wgp_mp = reshape((MPC(m).Wpg(p,:,:)) * MPC(m).M(p), [gridDimX gridDimY]);
        dxg_xp = Xg - reshape(MPC(m).X(p,:), [1,1,2]);
        BpDp_i = reshape(MPC(m).B(p,:,:), [2, 2]) * Dp_i;

        BpDp_i = reshape(BpDp_i, [1, 1, 2, 2]);

        temp = zeros(gridDimX, gridDimY, 2);
        temp(:,:,1) = BpDp_i(:,:,1,1).*dxg_xp(:,:,1) + BpDp_i(:,:,1,2).*dxg_xp(:,:,2);
        temp(:,:,2) = BpDp_i(:,:,2,1).*dxg_xp(:,:,1) + BpDp_i(:,:,2,2).*dxg_xp(:,:,2);
        MG(m).MOM = MG(m).MOM + wgp_mp .* (reshape(MPC(m).V(p,:), [1,1,2]) + temp);
    end

    %% 2. Compute grid velocities
    % non-vectorized, need to find vectorized solution
    for i=1:gridDimX
        for j=1:gridDimY
            if MG(m).M(i,j) == 0.0
                MG(m).V(i,j,:) = 0;
            else
                MG(m).V(i,j,:) = MG(m).MOM(i,j,:) / MG(m).M(i,j);
            end
        end
    end
    %% 3. Identify grid degrees of freedom
    [MG(m).degX, MG(m).degY] = find(MG(m).M);
    degX = MG(m).degX;
    degY = MG(m).degY;
    MG(m).GridDegrees = [degX, degY];
    numDegrees = size(degX);
    numDegrees = numDegrees(1);
    MG(m).numDegrees = numDegrees;
    MG(m).GInds1 = sub2ind([gridDimX gridDimY 2], degX, degY, ones(numDegrees, 1));
    MG(m).GInds2 = sub2ind([gridDimX gridDimY 2], degX, degY, 2*ones(numDegrees, 1));
    GInds1 = MG(m).GInds1;
    GInds2 = MG(m).GInds2;
    GridDegrees = MG(m).GridDegrees;
    %% 4. Compute explicit grid forces
    % compute determinant
    Jp = MPC(m).F(:,1,1).*MPC(m).F(:,2,2) - MPC(m).F(:,1,2).*MPC(m).F(:,2,1);

    % compute F inverse transpose
    Fp_it = NaN(MPC(m).N, 2, 2);
    for p=1:MPC(m).N
        Fp_it(p,:,:) = inv(reshape(MPC(m).F(p,:,:), [2 2]))'; 
    end

    % compute F transpose
    Fp_t = NaN(MPC(m).N, 2, 2);
    for p=1:MPC(m).N
        Fp_t(p,:,:) = reshape(MPC(m).F(p,:,:), [2 2])';
    end

    % compute Piola Kirchoff stress
    Pp = MPC(m).mew * (MPC(m).F - Fp_it) + MPC(m).lam * log(Jp).*Fp_it;

    MG(m).force = zeros(gridDimX, gridDimY, 2);
    for k=1:numDegrees 
        i = GridDegrees(k,1);
        j = GridDegrees(k,2);

        % element wise matrix*matrix multiplication
        temp1 = NaN(MPC(m).N,2,2);
        temp1(:,1,1) = Pp(:,1,1).*Fp_t(:,1,1) + Pp(:,1,2).*Fp_t(:,2,1);
        temp1(:,1,2) = Pp(:,1,1).*Fp_t(:,1,2) + Pp(:,1,2).*Fp_t(:,2,2);
        temp1(:,2,1) = Pp(:,2,1).*Fp_t(:,1,1) + Pp(:,2,2).*Fp_t(:,2,1);
        temp1(:,2,2) = Pp(:,2,1).*Fp_t(:,1,2) + Pp(:,2,2).*Fp_t(:,2,2);

        % element wise matrix * vector multiplication
        temp2 = NaN(MPC(m).N,2);
        temp2(:,1) = temp1(:,1,1).*MPC(m).Wpg_grad(:,i,j,1) + temp1(:,1,2).*MPC(m).Wpg_grad(:,i,j,2);
        temp2(:,2) = temp1(:,2,1).*MPC(m).Wpg_grad(:,i,j,1) + temp1(:,2,2).*MPC(m).Wpg_grad(:,i,j,2);

        MG(m).force(i,j,:) = -sum(MPC(m).VOL0.*temp2,1);
    end
    %% 5. Grid velocity update

    MG(m).V(GInds1) = MG(m).V(GInds1) + dt*MG(m).force(GInds1) ./ MG(m).M(GInds1);
    MG(m).V(GInds2) = MG(m).V(GInds2) + dt*MG(m).force(GInds2) ./ MG(m).M(GInds1) - dt*9.8;
end

% break here because we want to compute collisions for all material bodies
% at the same time

%% 5.2 Compute collisions between material point clouds
% find nodes where velocities are hitting eachother
[collisionX, collisionY] = find(sum(MG(1).V .* MG(2).V, 3) < 0.0);
colDegrees = size(collisionX);
colDegrees = colDegrees(1);
GColInds1 = sub2ind([gridDimX gridDimY 2], collisionX, collisionY, ones(colDegrees, 1));
GColInds2 = sub2ind([gridDimX gridDimY 2], collisionX, collisionY, 2*ones(colDegrees, 1));

% conservation of momentum: m1*v1i + m2*v2i = m1*v1f + m2*v2f
% v1f = m1
massDif = MG(1).M(GColInds1) - MG(2).M(GColInds1);
massSum = MG(1).M(GColInds1) + MG(2).M(GColInds1);
massTerm1 = massDif ./ massSum;
V1f.x = massTerm1 .* MG(1).V(GColInds1) + 2*MG(2).M(GColInds1) ./ massSum .*MG(2).V(GColInds1);
V1f.y = massTerm1 .* MG(1).V(GColInds2) + 2*MG(2).M(GColInds1) ./ massSum .*MG(2).V(GColInds2);

V2f.x = 2*MG(1).M(GColInds1) ./ massSum .*MG(1).V(GColInds1) - massTerm1 .* MG(2).V(GColInds1);
V2f.y = 2*MG(1).M(GColInds1) ./ massSum .*MG(1).V(GColInds2) - massTerm1 .* MG(2).V(GColInds2);

MG(1).V(GColInds1) = V1f.x;
MG(1).V(GColInds2) = V1f.y;
MG(2).V(GColInds1) = V2f.x;
MG(2).V(GColInds2) = V2f.y;

% for m=1:N_MPC
%     MG(m).V(GColInds1) = -MG(m).V(GColInds1);
%     MG(m).V(GColInds2) = -MG(m).V(GColInds2);
% end

% for m=1:N_MPC
%     numDegrees = MG(m).numDegrees;
%     GridDegrees = MG(m).GridDegrees;
%     GInds1 = MG(m).GInds1;
%     GInds2 = MG(m).GInds2;
%     degX = MG(m).degX;
%     degY = MG(m).degY;
%     %% 5.5 Compute collisions with grid
%     % compute grid node to border collisions
%     for k=1:numDegrees
%         i = GridDegrees(k,1);
%         j = GridDegrees(k,2);
% 
%         newNodePos = [i*h; j*h] + dt * reshape(MG(m).V(i,j,:), [2 1 1]);
% 
%         % left right borders
%         if (newNodePos(1) < BORDER_MIN_X) || (newNodePos(1) > BORDER_MAX_X)
%             MG(m).V(i, j, 1) = 0;
%             MG(m).V(i, j, 2) = MG(m).V(i, j, 2) * NONSTICKY;
%         end
% 
%         % top bottom borders
%         if (newNodePos(2) < BORDER_MIN_Y) || (newNodePos(2) > BORDER_MAX_Y)
%             MG(m).V(i, j, 1) = MG(m).V(i, j, 1) * NONSTICKY;
%             MG(m).V(i, j, 2) = 0;
%         end
%     end
% end

for m=1:N_MPC
    numDegrees = MG(m).numDegrees;
    GridDegrees = MG(m).GridDegrees;
    GInds1 = MG(m).GInds1;
    GInds2 = MG(m).GInds2;
    degX = MG(m).degX;
    degY = MG(m).degY;
    %% 6. Update particle deformation gradient
    for p=1:MPC(m).N
        PGInds1 = sub2ind([MPC(m).N gridDimX gridDimY 2], p*ones(numDegrees,1), degX, degY, 1*ones(numDegrees,1));
        PGInds2 = sub2ind([MPC(m).N gridDimX gridDimY 2], p*ones(numDegrees,1), degX, degY, 2*ones(numDegrees,1));

        temp = NaN(2,2);
        temp(1,1) = sum(MG(m).V(GInds1) .* MPC(m).Wpg_grad(PGInds1));
        temp(1,2) = sum(MG(m).V(GInds1) .* MPC(m).Wpg_grad(PGInds2));
        temp(2,1) = sum(MG(m).V(GInds2) .* MPC(m).Wpg_grad(PGInds1));
        temp(2,2) = sum(MG(m).V(GInds2) .* MPC(m).Wpg_grad(PGInds2));

          MPC(m).F(p,:,:) = (eye(2) + dt * temp) * reshape(MPC(m).F(p,:,:), [2 2 1]);
    end

    %% 7. Grid to particle transfer
    for p=1:MPC(m).N

        PGInds = sub2ind([MPC(m).N gridDimX gridDimY], p*ones(numDegrees,1), degX, degY);

        MPC(m).V(p,1) = sum(MPC(m).Wpg(PGInds).*MG(m).V(GInds1));
        MPC(m).V(p,2) = sum(MPC(m).Wpg(PGInds).*MG(m).V(GInds2));

        dXgp1 = Xg(GInds1) - MPC(m).X(p,1);
        dXgp2 = Xg(GInds2) - MPC(m).X(p,2);
        MPC(m).B(p,1,1) = sum(MPC(m).Wpg(PGInds).*MG(m).V(GInds1).*dXgp1);
        MPC(m).B(p,1,2) = sum(MPC(m).Wpg(PGInds).*MG(m).V(GInds1).*dXgp2);
        MPC(m).B(p,2,1) = sum(MPC(m).Wpg(PGInds).*MG(m).V(GInds2).*dXgp1);
        MPC(m).B(p,2,2) = sum(MPC(m).Wpg(PGInds).*MG(m).V(GInds2).*dXgp2);
    end
    %% 8. Particle advection
    MPC(m).X = MPC(m).X + dt*MPC(m).V;
end
frameTime = toc(frameStart);
fprintf("MPM timestep took %f seconds.\n", frameTime);

%% Plot/Animation
tic;

% draw the updated material point clouds
for m = 1:N_MPC
    delete(scatterPoints(m));
    scatterPoints(m) = scatter(MPC(m).X(:, 1), MPC(m).X(:,2), pointSize, 'filled', MPC(m).Color);
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