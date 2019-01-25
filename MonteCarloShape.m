clc;
clear all;


f = @(x,y) 2 - sqrt((x-5).^2 + (y-5).^2);

gridDimX = 20;
gridDimY = 20;
gridX0 = 1;
gridY0 = 1;

Np = 300;
p = 1;
Xp = NaN(Np, 2);
while p<=Np
    randX = gridX0 + rand() * (gridDimX - gridX0);
    randY = gridY0 + rand() * (gridDimY - gridY0);
    temp = f(randX, randY);
    if temp >= 0
        Xp(p,:) = [randX, randY];
        p = p + 1;
    end
end

%% Plot point cloud
pointSize = 3;
fig = figure;
plot([3 18 18 3 3], [3 3 18 18 3], 'k');
hold on;
scatter(Xp(:, 1), Xp(:,2), pointSize, 'red');
axis([1 20 1 20]);
xticks(1:20);
yticks(1:20);
grid on;
set(fig, 'Position', [50, 50, 700, 700]);