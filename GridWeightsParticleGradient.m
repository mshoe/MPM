function [W_grad] = GridWeightsParticleGradient(Xp, gridDimX, gridDimY, h)
%GRIDWEIGHTSPARTICLEGRADIENT Calculates the gradients of the weights wij
%for the given particle
%   Detailed explanation goes here
N_prime = @CubicBSplineDerivative;
N = @CubicBSpline;

xp = Xp(1);
yp = Xp(2);

W_grad = NaN(gridDimX, gridDimY, 2);
[IndX, IndY] = find(ones(gridDimX, gridDimY));
IndsX = sub2ind([gridDimX, gridDimY, 2], IndX, IndY, ones(size(IndX)));
IndsY = sub2ind([gridDimX, gridDimY, 2], IndX, IndY, 2*ones(size(IndX)));

W_grad(IndsX) = 1/h * N_prime(1/h*(xp-IndX)) .* N(1/h*(yp-IndY));
W_grad(IndsY) = 1/h * N(1/h*(xp-IndX)) .* N_prime(1/h*(yp-IndY));

% for i=1:gridDimX
%     for j=1:gridDimY
%         xij = i;
%         yij = j;
%         W_grad(i, j, 1) = 1/h*N_prime(1/h*(xp-xij)) * N(1/h*(yp-yij));
%         W_grad(i, j, 2) = N(1/h*(xp-xij)) * 1/h*N_prime(1/h*(yp-yij));
%     end
% end
end

