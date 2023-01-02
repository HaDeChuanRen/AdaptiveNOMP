clc; clear; close all;

rng(54);
c = 3e8;

% signal size parameter
Nx = 256;

% radar parameter
T_idle = 100e-6;
T_ramp = 60e-6;
T_circle = T_idle + T_ramp;
Fre_start = 77e9;
Slope_fre = 29.982e12;
Fs = 10e6;
Ts = 1 / Fs;
lambda_cw = c / Fre_start;
Rx_interval = lambda_cw / 2;

% the radial veclocity of target belongs to [-velocity_max, velocity_max]
velocity_max = c / (4 * Fre_start * T_circle);

% the radial range of target belongs to [0, range_max]
range_max = (c * Fs) / (2 * Slope_fre);
tau_max = 2 * range_max / c;
N_max = ceil(tau_max * Fs);

% the status of targets
K_targets = 4;
sigma_n = 1;
% SNR_all = 5 * rand(1, K_targets) + 1;
SNR_dB = 16;
% gain_alltargets = randn(K_targets, 1) + 1j * randn(K_targets, 1);
gain_allK = sqrt(10 .^ (SNR_dB / 10) * sqrt(sigma_n)) .* exp(1j * 2 * pi * rand(K_targets, 1));
% gain_allK = zeros(1, K_targets);

% set the variable to be estimated and get the corresponding angular frequency
velocity_set = 2 * velocity_max * rand(1, K_targets) -  velocity_max;
range_set = range_max * rand(1, K_targets);
theta_set = pi * rand(1, K_targets) - pi / 2;

% omega_x = 4 * pi * (Fre_start * velocity_set + Slope_fre * range_set) * Ts / c;
% omega_y = 4 * pi * Fre_start * velocity_set * T_circle / c;

omega_x = 1 : 4;
% omega_y = 5 : 8;

% omega_true = [omega_x; omega_y];
% omega_true = [1, 2, 3, 4; 5, 6, 7, 8];

array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);
y_vec = zeros(Nx, 1);

for k_idx = 1 : K_targets
    y_veck = array_Fun(omega_x(k_idx), Nx);
    y_vec = gain_allK(k_idx) * y_veck + y_vec;
end

y_vec = y_vec + sqrt(1 / 2) *(randn(Nx, 1) + 1j * randn(Nx, 1));

% gamma_mnomp = [4, 4];
K_known = K_targets;
% [omega_list, gain_list, ~] = MNOMP2D_Kknown_serial(y_matrix, gamma_mnomp, K_known);


guard_n = 2;
guard_m = 3;

training_n = 5;
training_m = 6;

guard_training_size = [guard_n, guard_m, training_n, training_m];

P_false_nominal = 1e-2;

gamma_cfar = [1, 1];
% gamma_mnomp = [4, 4, 4];
% MC = 1000;
% alpha_set = alpha_cal_2D_ca(gamma_cfar, P_false_nominal, Nx, My, guard_training_size, MC);
K_max = K_targets + 2;

% [omega_list, gain_list, y_residue_matrix, overthr_alldB] = MNOMP2D_CA_alpha_serial(y_matrix, gamma_cfar, gamma_mnomp, guard_training_size, alpha_set, K_max);

% p_fa_CFAR = 1e-2 / (Nx * My);
p_oe = 1e-2;
% N_r = (2 * training_n + 1) * (2 * training_m + 1) - (2 * guard_n + 1) * (2 * guard_m + 1);

N_r = 40;
S = eye(Nx);
alpha_hat = alpha_Poe(p_oe, Nx, N_r);
% [omegaList_CFAR, gainListCFAR, ~, Threshold_collect] =...
% MNOMP_CFAR_alpha(y_vec, S, alpha_hat, N_r, K_max);
[omegaList_CFAR, gainListCFAR, ~, Threshold_collect] = ...
MNOMP_CFARtr(y_vec, S, alpha_hat, N_r, K_max);

% S_snap = 1;
% tau = sigma_n * chi2inv((1 - p_oe) ^ (1 / Nx), 2 * S_snap) / 2;
% [omegaList_tau, gainList_tau, y_residue_matrix] =...
% MNOMP(y_vec, S, tau, K_max);



















