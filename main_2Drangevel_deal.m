clc; clear; close all;


% read the ADC data
orginal_path = 'C:\study\MNOMP_CFAR\4program\Matlab\data\20220506exp';
exp_type = '\01people1';
exp_serial = '\02';

filename = [orginal_path, exp_type, exp_serial, '\adc_data.bin'];
data_cube = readadc(filename);

% radar parameter
c = 3e8;
T_idle = 100e-6;
T_ramp = 60e-6;
T_circle = T_idle + T_ramp;
Fre_start = 77e9;
Slope_fre = 29.982e12;
Fs = 10e6;
Ts = 1 / Fs;
lambda_cw = c / Fre_start;
Rx_interval = lambda_cw / 2;


% set the observation scene
Nx = 64;
My = 32;
NM_num = Nx * My;

N_start = 1;
M_start = 17;
L_start = 1;
ymat = squeeze(data_cube(N_start : (N_start + Nx - 1), ...
M_start : (M_start + My - 1), L_start));

% algorithm parameter set
guard_n = 4;
guard_m = 2;
training_n = 6;
training_m = 6;
guard_training_size2D = [guard_n, guard_m, training_n, training_m];


K_max = 20;
P_oe2D = 1e-3;
N_r2D = (2 * training_n + 1) * (2 * training_m + 1) - ...
(2 * guard_n + 1) * (2 * guard_m + 1);

% NOMP method analysis
alpha_set2D = alpha_PoebyS(P_oe2D, N_r2D, NM_num);
[omegaList, gainList, ~, Threshold_collect] = ...
NOMP2D_CFAR(ymat, alpha_set2D, guard_training_size2D, K_max);


% NOMP2D_CFAR(y_matrix, p_fa_CFAR, guard_training_size2D, K_max);
K_est = length(gainList);

omegax_hat = omegaList(:, 1);
omegay_hat = wrapToPi(omegaList(:, 2));
range_hat = (c * omegax_hat) / (4 * pi * Ts * Slope_fre);
velocity_hat = (c * omegay_hat) / (4 * pi * Fre_start * T_circle);
target_amp = 20 * log10(abs(gainList));
target_threshold = 10 * log10(abs(Threshold_collect));

% M_velocity = 1024;
% velocity_idx_long = linspace(- velocity_max, velocity_max, M_velocity);

% target_amp = nan(N_range, M_velocity);
% target_threshold = nan(N_range, M_velocity);

% target_idxn =  round(N_range * omegax_hat / (2 * pi));
% target_idxm =  round(M_velocity * (omegay_hat + pi) / (2 * pi));

% for k_idx = 1 : K_est
%     target_amp(target_idxn(k_idx) + 1, target_idxm(k_idx)) = 20 * log10(abs(gainList(k_idx)));
%     target_threshold(target_idxn(k_idx) + 1, target_idxm(k_idx)) = 10 * log10(abs(Threshold_collect(k_idx)));
% end

% [range_meshidx, velocity_meshidx] = meshgrid(range_idx_long, velocity_idx_long);

figure;
stem3(range_hat, velocity_hat, target_amp');
hold on;
stem3(range_hat, velocity_hat, target_threshold', 'r*');
xlim([0 25])
xlabel('range(m)');
ylabel('velocity(m/s)')
zlabel('amplitude(dB)')
legend('targets amplitude', 'thereshold')






