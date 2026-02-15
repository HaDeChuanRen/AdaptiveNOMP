clc; clear; close all;

Nx = 512;
n_var = 0 : Nx;
alpha_set = 12;
N_r = 60;

resvec_PFAln = zeros(Nx + 1, 1);
PFAvec_res = zeros(Nx + 1, 1);

for n_idx = 0 : Nx
    resvec_PFAln(n_idx + 1) = sum(log(1 : Nx)) - sum(log(1 : n_idx)) - ...
        sum(log(1 : Nx - n_idx)) - N_r * log(1 + n_idx * alpha_set / N_r);
end
Pfa_N = sum((- 1) .^ n_var' .* exp(resvec_PFAln));

fun_Poe = @(xvar) exp(Nx * log(1 - exp(- alpha_set / (2 * N_r) * xvar))+ ...
    (N_r - 1) * log(xvar) - xvar / 2 - N_r * log(2) - sum(log(1 : N_r - 1)));
Pfa_int = 1 - integral(fun_Poe, 0, Inf);

lw = 1.2;
fsz = 12;
msz = 8;

figure;
stem(n_var, resvec_PFAln, 'Markersize', msz)
xlim([0 512])
xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('$G_n$', 'Interpreter', 'latex', 'Fontsize', fsz)



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

% sigma_n = 1;
% 
% sigma_vec = sigma_n * exprnd(1, Nx);
% y_noise = sqrt(sigma_vec / 2) .* (randn(Nx, 1) + 1j * randn(Nx, 1));
% 
% norm(y_noise) / Nx





