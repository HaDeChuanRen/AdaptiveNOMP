clc; clear; close all;

Nx = 256000;
mu_sigma = 1;
sigma_n = 1;
sigma_vec = sigma_n * exprnd(mu_sigma, [Nx, 1]);
y_noise = sqrt(sigma_vec / 2) .* (randn(Nx, 1) + 1j * randn(Nx, 1));
% mean(sigma_vec)


