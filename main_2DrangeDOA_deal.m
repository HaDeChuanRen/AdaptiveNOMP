clc; clear; close all;


% read the ADC data
orginal_path = 'D:\XuMenghuai\FMCW mmwave range and Doppler estimation 202106\4program\Matlab\data\20220506exp\20220506exp';
exp_type = '\02people2';
exp_serial = '\03';

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
Nx = 256;
Lz = 4;
NL_num = Nx * Lz;

N_start = 1;
M_start = 17;
L_start = 1;
ymat = squeeze(data_cube(N_start : (N_start + Nx - 1), M_start, :));

% algorithm parameter set
guard_n = 4;
guard_l = 0;
training_n = 10;
training_l = 1;
guard_training_size2D = [guard_n, guard_l, training_n, training_l];


K_max = 20;
P_oe2D = 1e-3;
N_r2D = (2 * training_n + 1) * (2 * training_l + 1) - ...
(2 * guard_n + 1) * (2 * guard_l + 1);
alpha_set2D = alpha_PoebyS(P_oe2D, N_r2D, NL_num);

% NOMP method analysis

[omega_DOA, gain_DOA, ~, Threshold_DOA] = ...
NOMP2D_CFAR(ymat, alpha_set2D, guard_training_size2D, K_max);

omegax_hat = omega_DOA(:, 1);
omegaz_hat = wrapToPi(omega_DOA(:, 2));
range_hat = (c * omegax_hat) / (4 * pi * Ts * Slope_fre);
theta_hat = asin((c * omegaz_hat) / (2 * pi * Fre_start * Rx_interval));

xloc_hat = range_hat .* sin(theta_hat);
yloc_hat = range_hat .* cos(theta_hat);

figure;
stem3(xloc_hat, yloc_hat, 20 * log10(abs(gain_DOA)));
hold on;
plot3(xloc_hat, yloc_hat, 10 * log10(abs(Threshold_DOA)), 'r*');
xlim([-5 5])
ylim([0 6])
xlabel('xaxis(m)');
ylabel('yaxis(m)');
zlabel('amplitude(dB)')
legend('targets amplitude', 'thereshold')
title('location estimation')




