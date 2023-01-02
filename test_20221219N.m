clc; clear; close all;

Nx = 256;
% Nr = 50;
% xvar = 0 : 1000;
% alpha_set = 10.22;

% fun_ln = N * log(1 - exp(- alpha_set / (2 * Nr) * xvar)) + ...
%     (Nr - 1) * log(xvar) - xvar / 2 - Nr * log(2) - sum(log(1 : Nr - 1));

% % plot the result
% lw = 2;
% fsz = 12;
% msz = 8;

% figure;
% plot(xvar, fun_ln, 'r', 'Linewidth', lw, 'Markersize', msz);
% xlabel('x', 'Fontsize', fsz)
% ylabel('f(x)', 'Fontsize', fsz)

sigma_n = 1;

sigma_vec = sigma_n * exprnd(1, Nx);
y_noise = sqrt(sigma_vec / 2) .* (randn(Nx, 1) + 1j * randn(Nx, 1));

norm(y_noise) / Nx





