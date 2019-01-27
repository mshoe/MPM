clc;
clear all;

gridDimX = 20;
gridDimY = 20;
h = 1;
d = h/4;
gridX0 = 1;
gridY0 = 1;

[Xp, Np] = PointCloud(d, gridDimX, gridDimY, gridX0, gridY0);

% mass of material points
pointMass = 1;
Mp = pointMass * ones(Np, 1);

% intial weights
Wpg = GridWeights(Xp, Np, gridDimX, gridDimY, h); 
% initial volumes
VOL0p = InitialParticleVolumes(Wpg, Mp, Np, gridDimX, gridDimY, h);

%% Plot point cloud
pointSize = 4/d;
fig = figure;%('visible','off');

plot([3 18 18 3 3 18], [3 3 18 18 3 3], 'w', 'LineWidth',7);
set(gca, 'Color', 'k');
ax = gca;
ax.GridColor = [1.0, 1.0, 1.0];
hold on;
props = {'LineStyle','none','Marker','o','MarkerEdge','b','MarkerSize',6};
MPs = scatter(Xp(:, 1), Xp(:,2), pointSize, 'filled', 'red');
%MPs = line([X0p(:,1), X0p(:,1)], [X0p(:,2), X0p(:,2)], props{:});
axis([1 20 1 20]);
xticks(1:20);
yticks(1:20);
grid on;
set(fig, 'Position', [50, 50, 700, 700]);