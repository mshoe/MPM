function [W_grad] = GridWeightsParticleGradient(Xp, gridDimX, gridDimY, gridX0, gridY0, h)
%GRIDWEIGHTSPARTICLEGRADIENT Calculates the gradients of the weights wij
%for the given particle
%   Detailed explanation goes here
N_prime = @CubicBSplineDerivative;
N = @CubicBSpline;

xp = Xp(1);
yp = Xp(2);

W_grad = NaN(gridDimX, gridDimY, 2);

for i=gridX0:gridDimX
    for j=gridY0:gridDimY
        xij = gridX0 + i - 1;
        yij = gridY0 + j - 1;
        W_grad(i, j, 1) = 1/h*N_prime(1/h*(xp-xij)) * N(1/h*(yp-yij));
        W_grad(i, j, 2) = N(1/h*(xp-xij)) * 1/h*N_prime(1/h*(yp-yij));
    end
end
end

