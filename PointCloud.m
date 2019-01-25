function [Xp] = PointCloud(Np, gridDimX, gridDimY, gridX0, gridY0)
%POINTCLOUD Summary of this function goes here
%   Detailed explanation goes here


circle = @(p, c, r) sqrt((p.x-c.x).^2 + (p.y-c.y).^2) - r;
donut = @(p, c, r1, r2) max(-circle(p, c, r1), circle(p, c, r2));


c.x = 12;
c.y = 8;
r = 2;

i = 1;
Xp = NaN(Np, 2);
while i<=Np
    p.x = gridX0 + rand() * (gridDimX - gridX0);
    p.y = gridY0 + rand() * (gridDimY - gridY0);
    
    temp = circle(p, c, r);
    %temp = donut(p, c, r, r*2);
    if temp < 0
        Xp(i,:) = [p.x, p.y];
        i = i + 1;
    end
end
end

