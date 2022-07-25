clc; clear; close all;

orginal_path = 'C:\study\MNOMP_CFAR\4program\Matlab\data\20220506exp';
exp_type = '\01people1';
exp_serial = '\02';

filename = [orginal_path, exp_type, exp_serial, '\adc_data.bin'];

data_cube = readadc(filename);

K_max = 10;

Nx = 256;
My = 16;
Lz = 4;
size_D = [Nx, My, Lz];
S = eye(Nx);

NM_num = Nx * My;

y_tensor = data_cube(1 : Nx, 1 : My, :);

% P_FA_nominal = 1e-1;
p_fa_CFAR = 1e-4 / Nx;
MC = 1000;
gamma_cfar = [1, 1, 1];
gamma_cfar1D = 1;
gamma_mnomp = [4, 4, 4];
gamma_mnomp1D = 4;
guard_n = 4;
guard_m = 2;
guard_l = 1;

training_n = 30;
training_m = 3;
training_l = 2;
guard_training_size1D = [guard_n, training_n];
guard_training_size2D = [guard_n, guard_m, training_n, training_m];

P_oe = 0.01;
N_r = 60;
alpha_set = alpha_PoebyS(P_oe, Nx, N_r);

% alpha_set = alpha_cal_ca(gamma_cfar1D, P_FA_nominal, Nx, guard_training_size1D, MC);

y_vector = data_cube(:, 128, 1);
y_matrix = squeeze(data_cube(1 : 128, 65 : (65 + My), 1));


% [omega_list, gain_list, y_residue_vector] = MNOMP_Kaware_serial(y_tensor, gamma_mnomp, K_max);
% [omega_list, gain_list, y_residue, overthr_alldB, sigma_hat] = MNOMP_CA_alpha_serial(y_vector, gamma_cfar1D, gamma_mnomp, guard_training_size1D, alpha_set, K_max);
% [omega_list, gain_list, y_residue, overthr_alldB, sigma_hat] = MNOMP_CA_backward_serial(y_vector, gamma_cfar1D, gamma_mnomp1D, guard_training_size1D, alpha_set, K_max);
% [omega_list, gain_list, residueList_cfar] = MNOMP_CFAR(y_vector, S, alpha_set, guard_training_size1D, K_max);
% [omega_list, gain_list, residueList_cfar] = MNOMP(y_vector, S, tau);
% [omega_list, gain_list, residueList, Threshold_collect] = NOMP1D_CFAR(y_vector, p_fa_CFAR, guard_training_size1D, K_max);

[omega_list, gain_list, residueList, Threshold_collect] = ...
MNOMP_CFAR_alpha(y_vector, S, alpha_set, N_r, K_max);


omegax_hat = omega_list;
% omegay_hat = wrapToPi(omega_list(2, :));
% omegaz_hat = wrapToPi(omega_list(3, :));

c = 3e8;
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

range_max = (c * 2 * pi) / (4 * pi * Ts * Slope_fre);
velocity_max = (c * pi) / (4 * pi * Fre_start * T_circle);

range_idx = linspace(0, range_max, Nx);
figure;
subplot(2, 1, 1)
plot(range_idx, 20 * log10(abs(fft(squeeze(y_vector)) / sqrt(Nx))));
xlabel('range(m)')
ylabel('amplitude(dB)')
title('FFT result of IF')
xlim([0, range_max])
% (c * 2 * pi) / (4 * pi * Ts * Slope_fre)
range_hat = (c * omegax_hat) / (4 * pi * Ts * Slope_fre)
% velocity_hat = (c * omegay_hat) / (4 * pi * Fre_start * T_circle);
% theta_hat = asin((c * omegaz_hat) / (2 * pi * Fre_start * Rx_interval));

N_range = 2048;
range_idx_long = linspace(0, range_max, N_range);
target_amp = nan(N_range, 1);
target_threshold = nan(N_range, 1);


target_idx = round(N_range * omega_list / (2 * pi));
target_amp(target_idx + 1) = 20 * log10(abs(gain_list));
target_threshold(target_idx + 1) = 10 * log10(abs(Threshold_collect));


subplot(2, 1, 2)
stem(range_idx_long, target_amp);
hold on;
plot(range_idx_long, target_threshold, 'r*');
legend('targets amplitude', 'thereshold')
xlabel('range(m)')
ylabel('amplitude(dB)')
title('NOMP-CFAR result of IF')
xlim([0, range_max])


N_r2D = (2 * training_n + 1) * (2 * training_m + 1) - ...
(2 * guard_n + 1) * (2 * guard_m + 1);
alpha_set2D = alpha_PoebyS(P_oe, N_r2D, NM_num);
[omegaList, gainList, ~, Threshold_collect] = ...
NOMP2D_CFAR(y_matrix, alpha_set2D, guard_training_size2D, K_max);


% NOMP2D_CFAR(y_matrix, p_fa_CFAR, guard_training_size2D, K_max);
K_est = length(gainList);

omegax_hat = omegaList(:, 1);
omegay_hat = wrapToPi(omegaList(:, 2));
range_hat = (c * omegax_hat) / (4 * pi * Ts * Slope_fre);
velocity_hat = (c * omegay_hat) / (4 * pi * Fre_start * T_circle);

M_velocity = 1024;
velocity_idx_long = linspace(- velocity_max, velocity_max, M_velocity);

target_amp = nan(N_range, M_velocity);
target_threshold = nan(N_range, M_velocity);

target_idxn =  round(N_range * omegax_hat / (2 * pi));
target_idxm =  round(M_velocity * (omegay_hat + pi) / (2 * pi));

for k_idx = 1 : K_est
    target_amp(target_idxn(k_idx) + 1, target_idxm(k_idx)) = 20 * log10(abs(gainList(k_idx)));
    target_threshold(target_idxn(k_idx) + 1, target_idxm(k_idx)) = 10 * log10(abs(Threshold_collect(k_idx)));
end

[range_meshidx, velocity_meshidx] = meshgrid(range_idx_long, velocity_idx_long);

figure;
stem3(range_meshidx, velocity_meshidx, target_amp');
hold on;
stem3(range_meshidx, velocity_meshidx, target_threshold', 'r*');
xlim([0 25])
xlabel('range(m)');
ylabel('velocity(m/s)')
zlabel('amplitude(dB)')
legend('targets amplitude', 'thereshold')

