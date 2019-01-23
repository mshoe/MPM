clear all;
clc;

%% Grid initialization

% bounds
% x,y bounds
bounds = [3 18];

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
% number of material points
Np = 300;

% reference configuration of material points
%X0p = zeros(Np, 2);%, 'gpuArray');
BoxBounds = [5 8];
X0p = BoxBounds(1) + rand(Np, 2) * (BoxBounds(2) - BoxBounds(1));

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

%% MPM final initializations

% intial weights
Wpg = GridWeights(X0p, Np, gridDimX, gridDimY, gridX0, gridY0, h); 

% initial volumes
VOL0p = InitialParticleVolumes(Wpg, Mp, Np, gridDimX, gridDimY, h);

% initial APIC B matrices
Bp = zeros(Np, 2, 2);

% this matrix doesn't change? APIC D
Dp_i = inv(1/3 * h^2 * eye(2));

%% MPM Algorithm
% ************************************************************** %
timeSteps = 1000;
dt = 1/60;

pointSize = 3;
fig = figure;
plot([3 18 18 3 3], [3 3 18 18 3], 'k');
hold on;
scatter(X0p(:, 1), X0p(:,2), pointSize, 'red');
axis([1 20 1 20]);
xticks(1:20);
yticks(1:20);
grid on;
set(fig, 'Position', [50, 50, 700, 700]);
M(1) = getframe;

v = VideoWriter('mpm_test.avi');
open(v);

for t=1:timeSteps
%% 1. Particle to grid transfer (using APIC)
tic;
% first get the grid weights for each material point
Wpg = GridWeights(Xp, Np, gridDimX, gridDimY, gridX0, gridY0, h); 
Wpg_grad = GridWeightsGradient(Xp, Np, gridDimX, gridDimY, gridX0, gridY0, h);

% then transfer the masses of the mps to the grid
Mg = zeros(gridDimX, gridDimY);
for p=1:Np
    Mg = Mg + reshape(Wpg(p,:,:), [gridDimX, gridDimY]) * Mp(p);
end

% transfer the momentum of the mps to the grid using APIC formula
MOMg = zeros(gridDimX, gridDimY, 2);
for p=1:Np
    
    wgp_mp = squeeze(Wpg(p,:,:)) * Mp(p);
    dxg_xp = Xg - reshape(Xp(p,:), [1,1,2]);
    BpDp_i = reshape(Bp(p,:,:), [2, 2]) * Dp_i;
    
    BpDp_i = reshape(BpDp_i, [1, 1, 2, 2]);
    
    % super slow non-vectorized code, need to find vectorized solution
    temp = zeros(gridDimX, gridDimY, 2);
    temp(:,:,1) = BpDp_i(:,:,1,1).*dxg_xp(:,:,1) + BpDp_i(:,:,1,2).*dxg_xp(:,:,2);
    temp(:,:,2) = BpDp_i(:,:,2,1).*dxg_xp(:,:,1) + BpDp_i(:,:,2,2).*dxg_xp(:,:,2);
%     for i=1:gridDimX
%         for j = 1:gridDimY
%             temp(i, j, :) = BpDp_i * squeeze(dxg_xp(i, j, :));
%         end
%     end
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

squeeze(Vg(6,7,:))

%% 3. Identify grid degrees of freedom

GridDegrees = [];
numDegrees = 0;
for i=1:gridDimX
    for j=1:gridDimY
        if Mg(i, j) ~= 0.0
            GridDegrees = cat(1, GridDegrees, [i j]);
            numDegrees = numDegrees + 1;
        end
    end
end

%% 4. Compute explicit grid forces

% compute determinant
Jp = Fp(:,1,1).*Fp(:,2,2) - Fp(:,1,2).*Fp(:,2,1);

% compute F inverse transpose
Fp_it = NaN(Np, 2, 2);
for p=1:Np
    Fp_it(p,:,:) = inv(squeeze(Fp(p,:,:)))'; 
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
%     temp1 = NaN(Np,2,2);
%     temp1(:,1,1) = Pp(:,1,1).*Fp_t(:,1,1) + Pp(:,1,2).*Fp_t(:,2,1);
%     temp1(:,1,2) = Pp(:,1,1).*Fp_t(:,1,2) + Pp(:,1,2).*Fp_t(:,2,2);
%     temp1(:,2,1) = Pp(:,2,1).*Fp_t(:,1,1) + Pp(:,2,2).*Fp_t(:,2,1);
%     temp1(:,2,2) = Pp(:,2,1).*Fp_t(:,1,2) + Pp(:,2,2).*Fp_t(:,2,2);
%     
%     % element wise matrix * vector multiplication
%     temp2 = NaN(Np,2);
%     temp2(:,1) = temp1(:,1,1).*Wpg_grad(:,1) + temp1(:,1,2).*Wpg_grad(:,2);
%     temp2(:,2) = temp1(:,2,1).*Wpg_grad(:,1) + temp1(:,2,2).*Wpg_grad(:,2);
%     
%     force_g(i,j,:) = -sum(temp2,1);
    for p=1:Np
        temp = squeeze(Pp(p,:,:)) * squeeze(Fp_t(p,:,:)) * squeeze(Wpg_grad(p,i,j,:));
        force_g(i,j,:) = force_g(i,j,:) - reshape(VOL0p(p)*temp, [1 1 2]);
    end
end
squeeze(force_g(7,6,:))
%% 5. Grid velocity update

for k=1:numDegrees
    i = GridDegrees(k,1);
    j = GridDegrees(k,2);
    Vg(i,j,:) = Vg(i,j,:) + dt*force_g(i,j,:) / Mg(i,j);
    Vg(i,j,2) = Vg(i,j,2) - dt * 9.8;
end

% compute grid node to border collisions
for k=1:numDegrees
    i = GridDegrees(k,1);
    j = GridDegrees(k,2);
    
    newNodePos = [i*h; j*h] + dt * squeeze(Vg(i,j,:));
    
    % left right borders
    if (newNodePos(1) < BORDER_MIN) || (newNodePos(1) > BORDER_MAX)
        Vg(i, j, 1) = 0;
    end
    
    % top bottom borders
    if (newNodePos(2) < BORDER_MIN) || (newNodePos(2) > BORDER_MAX)
        Vg(i, j, 2) = 0;
    end
end

%% 6. Update particle deformation gradient

for p=1:Np
    temp = zeros(2,2);
    for k=1:numDegrees
        i = GridDegrees(k,1);
        j = GridDegrees(k,2);

        temp = temp + squeeze(Vg(i,j,:)) * squeeze(Wpg_grad(p,i,j,:))';
    end
    Fp(p,:,:) = (eye(2) + dt * temp) * squeeze(Fp(p,:,:));
end

%% 7. Grid to particle transfer

for p=1:Np
    temp = [0 0];
    for k=1:numDegrees
        i = GridDegrees(k,1);
        j = GridDegrees(k,2);
        
        temp = temp + Wpg(p,i,j) * squeeze(Vg(i,j,:))';
    end
    Vp(p,:) = temp;
    %Vp(p,2) = Vp(p,2) - dt*9.8;
end

for p=1:Np
    temp = zeros(2);
    for k=1:numDegrees
        i = GridDegrees(k,1);
        j = GridDegrees(k,2);
        
        temp = temp + Wpg(p,i,j) * squeeze(Vg(i,j,:)) * squeeze(Xg(i,j,:) - reshape(Xp(p,:), [1 1 2]))';
    end
    Bp(p,:,:) = temp;
end

%% 8. Particle advection


Xp = Xp + dt*Vp;

%% Plot/Animation
t
toc;
clf;

plot([3 18 18 3 3], [3 3 18 18 3], 'k');
hold on;
scatter(Xp(:, 1), Xp(:,2), pointSize, 'red');
axis([1 20 1 20]);
xticks(1:20);
yticks(1:20);
grid on;
set(fig, 'Position', [50, 50, 700, 700]);
M(t+1) = getframe;

end

writeVideo(v, M);
close(v);