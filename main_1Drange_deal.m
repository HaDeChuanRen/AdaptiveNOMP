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
Smat_com = eye(Nx);

M_start = 17;
L_start = 1;
yvec = data_cube(:, M_start, L_start);

% algorithm parameter set
K_max = 10;
P_oe = 1e-3;
guard_n = 4;
training_n = 30;
N_r = training_n * 2;
alpha_set = alpha_PoebyS(P_oe, Nx, N_r);
% 10 * log10(alpha_set)

% CA-NOMP method analysis
tic;
[omega_list, gain_list, residueList, Threshold_collect] = ...
MNOMP_CFAR_alpha(yvec, Smat_com, alpha_set, N_r, K_max);
toc;


% CA-CFAR and FFT analysis
tic;
y_fftabs_vector = abs(fft(squeeze(yvec)) / sqrt(Nx));
prob_ind_ext = repmat(y_fftabs_vector .^ 2, [3, 1]);
cfar_detector = phased.CFARDetector('NumTrainingCells', 2 * training_n, 'NumGuardCells', 2 * guard_n);
cfar_detector.Method = 'CA';
% if strcmp(CFAR_method, 'OS')
%     N_r = (2 * training_n + 1) - (2 * guard_n + 1);
%     cfar_detector.Rank = round(N_r / 2);
% end
cfar_detector.ThresholdFactor = 'Auto';
cfar_detector.ProbabilityFalseAlarm = P_oe / Nx;
cfar_detector.ThresholdOutputPort = true;
cfar_detector.NoisePowerOutputPort = true;
[peak_grid, Threshold_CUT, sigma_hat] = cfar_detector(prob_ind_ext, (Nx + 1 : 2 * Nx)');
toc;

% figure;
% plot(peak_grid)
% T_judgement = res_inf_normSq_rot / Threshold_CUT - 1;

range_max = (c * 2 * pi) / (4 * pi * Ts * Slope_fre);
% velocity_max = (c * pi) / (4 * pi * Fre_start * T_circle);
range_idx = linspace(0, range_max, Nx);
range_hatfft = range_idx(peak_grid);
figure;
% subplot(2, 1, 1)
plot(range_idx, 20 * log10(y_fftabs_vector));
hold on;
stem(range_hatfft, 20 * log10(y_fftabs_vector(peak_grid)), 'bo');
plot(range_hatfft, 10 * log10(Threshold_CUT(peak_grid)), 'r*')
legend('FFT result', 'target', 'thereshold')
xlabel('range(m)')
ylabel('amplitude(dB)')
title('FFT result of IF')
xlim([0, range_max / 5])
% (c * 2 * pi) / (4 * pi * Ts * Slope_fre)

omegax_hat = omega_list;
range_hat = (c * omegax_hat) / (4 * pi * Ts * Slope_fre)
% subplot(2, 1, 2)
figure;
stem(range_hat, 20 * log10(abs(gain_list)));
hold on;
plot(range_hat, 10 * log10(abs(Threshold_collect)), 'r*');
legend('target', 'thereshold')
xlabel('range(m)')
ylabel('amplitude(dB)')
title('NOMP-CFAR result of IF')
xlim([0, range_max / 5])
























% figure;
% subplot(2, 1, 1)
% plot(range_idx, 20 * log10(abs(fft(squeeze(y_vector)) / sqrt(Nx))));
% xlabel('range(m)')
% ylabel('amplitude(dB)')
% title('FFT result of IF')
% xlim([0, range_max])
% % (c * 2 * pi) / (4 * pi * Ts * Slope_fre)
% range_hat = (c * omegax_hat) / (4 * pi * Ts * Slope_fre)
% % velocity_hat = (c * omegay_hat) / (4 * pi * Fre_start * T_circle);
% % theta_hat = asin((c * omegaz_hat) / (2 * pi * Fre_start * Rx_interval));

% N_range = 2048;
% range_idx_long = linspace(0, range_max, N_range);
% target_amp = nan(N_range, 1);
% target_threshold = nan(N_range, 1);


% target_idx = round(N_range * omega_list / (2 * pi));
% target_amp(target_idx + 1) = 20 * log10(abs(gain_list));
% target_threshold(target_idx + 1) = 10 * log10(abs(Threshold_collect));


% subplot(2, 1, 2)
% stem(range_idx_long, target_amp);
% hold on;
% plot(range_idx_long, target_threshold, 'r*');
% legend('targets amplitude', 'thereshold')
% xlabel('range(m)')
% ylabel('amplitude(dB)')
% title('NOMP-CFAR result of IF')
% xlim([0, range_max])




