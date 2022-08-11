clc; clear; close all;

rng(5);
MC = 5000;

% Define Scenario
Nx = 2048; % Length of Sinusoid
K = 8;
sigma_n = 1;              % noise variance sigma^2, instead of sigma
S_snap = 1;
SNR = 22;

M = Nx;
Smat_com = eye(M);

% algorithm parameters
K_max = 10;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;
N_r = 60;
alpha_set = alpha_PoebyS(P_oe, Nx, N_r, S_snap);
tau_set = sigma_n * chi2inv((1 - P_oe) ^ (1 / Nx), 2 * S_snap) / 2;

% statistical variable initialize
Equalmat_tau = zeros(MC, 1);
Equalmat_CA = zeros(MC, 1);
Equalmat_prob = zeros(MC, 1);

ymat_colloct = zeros(Nx, MC);

for mc_idx = 1 : MC
    [y, omega_true, gain_true] = create_yvector(K, S_snap, SNR, sigma_n, Smat_com);
    ymat_colloct(:, mc_idx) = y;
end

tic;
for mc_idx = 1 : MC
    [omegavec_tau, gainvec_tau, ~] = MNOMP(ymat_colloct(:, mc_idx), Smat_com, tau_set);
    if length(omegavec_tau) == K
        Equalmat_tau(mc_idx) = 1;
    end
end
time_NOMP = toc;


tic;
for mc_idx = 1 : MC
    % handle_waitbar = waitbar(mc_idx / MC);
    % [y, omega_true, gain_true] = create_yvector(K, S_snap, SNR, sigma_n, Smat_com);
    [omegavec_CA, gainvec_CA, ~] = ...
    MNOMP_CFAR_alpha(ymat_colloct(:, mc_idx), Smat_com, alpha_set, N_r, K_max);
    if length(omegavec_CA) == K
        Equalmat_CA(mc_idx) = 1;
    end
end
time_CANOMP = toc;

% tic;
% for mc_idx = 1 : MC
%     % handle_waitbar = waitbar(mc_idx / MC);
%     % [y, omega_true, gain_true] = create_yvector(K, S_snap, SNR, sigma_n, Smat_com);
%     [omegavec_prob, gainvec_prob, ~] = ...
%     MNOMP_alpha_prob(ymat_colloct(:, mc_idx), Smat_com, alpha_set, N_r, K_max);
%     if length(omegavec_prob) == K
%         Equalmat_prob(mc_idx) = 1;
%     end
% end
% time_probNOMP = toc;
% 
% time_probNOMP - time_CANOMP

% mean(Equalmat_CA)
% mean(Equalmat_prob)

