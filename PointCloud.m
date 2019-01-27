function [Xp, Np] = PointCloud(shape, d, gridDimX, gridDimY, gridX0, gridY0)
%POINTCLOUD
%   S (shape) is a struct with possible parameters:
%       string TYPE
%       c1 (c.x, c.y) center of shape
%       r1
%       r2 (r2 > r1)

%% Primitives
circle = @(p, c, r) sqrt((p.x-c.x).^2 + (p.y-c.y).^2) - r;
donut = @(p, c, r1, r2) max(-circle(p, c, r1), circle(p, c, r2));

union = @(sd1, sd2) min(sd1, sd2);

%% Definied demo shapes
demoCircle = @(s, p) circle(p, s.c, s.r);
demoDonut = @(s, p) donut(p, s.c, s.r1, s.r2);
twoCirclesUnion = @(s, p) union(circle(p,s.c1,s.r),circle(p,s.c2,s.r));

sdf = demoCircle;

if shape.type=="twoCirclesUnion"
    sdf = twoCirclesUnion;
elseif shape.type=="demoDonut"
    sdf = demoDonut;
end


Xp = [];
for i=gridX0:d:gridDimX
    for j=gridY0:d:gridDimY
        p.x = i;
        p.y = j;
        
        %temp = donut(p, c, r, r^2);
        temp = sdf(shape, p);
        if temp < 0
            Xp = [Xp; p.x p.y];
        end
    end
end
Np = size(Xp);
Np = Np(1);


% while i<=Np
%     p.x = gridX0 + rand() * (gridDimX - gridX0);
%     p.y = gridY0 + rand() * (gridDimY - gridY0);
%     
%     temp = circle(p, c, r);
%     %temp = donut(p, c, r, r*2);
%     if temp < 0
%         Xp(i,:) = [p.x, p.y];
%         i = i + 1;
%     end
% end
end

