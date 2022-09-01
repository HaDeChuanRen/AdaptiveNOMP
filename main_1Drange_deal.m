clc; clear; close all;


% read the ADC data
orginal_path = 'C:\study\MNOMP_CFAR\4program\Matlab\data\20220506exp';
exp_type = '\02people2';
exp_serial = '\01';

range_true = [2.95 + 0.14; 4.78 + 0.14];
K_true = length(range_true);
amplitude_idx = 60;

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
K_max = 32;
P_oe = 1e-2;
guard_n = 4;
training_n = 30;
N_r = training_n * 2;
alpha_set = alpha_PoebyS(P_oe, Nx, N_r);
sigma_set = 10 ^ (24 / 10);
tau_set = sigma_set * chi2inv((1 - P_oe) ^ (1 / Nx), 2) / 2;
% 10 * log10(alpha_set)

% CA-NOMP method analysis
tic;
[omega_list, gain_list, residueList, Threshold_collect] = ...
MNOMP_CFAR_alpha(yvec, Smat_com, alpha_set, N_r, K_max);
time_NOMPCFAR = toc;

tic;
[omegalist_tau, gainlist_tau, ~] = MNOMP(yvec, Smat_com, tau_set);
time_tau = toc;

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
time_CFAR = toc;

% Threshold_CUT ./ sigma_hat
% figure;
% plot(peak_grid)
% T_judgement = res_inf_normSq_rot / Threshold_CUT - 1;

range_max = (c * 2 * pi) / (4 * pi * Ts * Slope_fre);
% velocity_max = (c * pi) / (4 * pi * Fre_start * T_circle);
range_idx = linspace(0, range_max, Nx);
range_hatfft = range_idx(peak_grid);

lw = 2;
fsz = 12;
msz = 8;

bias = 10 * log10(Nx);

% figure;
% plot(20 * log10(y_fftabs_vector))

figure;
% subplot(2, 1, 1)
plot(range_idx, 20 * log10(y_fftabs_vector), '-k', 'Linewidth', lw);
hold on;
plot(range_idx, 10 * log10(Threshold_CUT), '-.r', 'Linewidth', lw)

% plot(range_true(1)*ones(size(amplitude_idx)), amplitude_idx,':r', range_true(2)*ones(size(amplitude_idx)), amplitude_idx,':r', 'Linewidth', lw);
% stem(range_true(1), amplitude_idx * ones(K_true, 1), ':.r', 'Linewidth', lw,'DisplayName','True', 'Fontsize', fsz)
plot(range_hatfft, 20 * log10(y_fftabs_vector(peak_grid)) , 'bo', 'Linewidth', lw, 'Markersize', msz);
stem(range_true, amplitude_idx * ones(K_true, 1) , ':.m', 'Linewidth', lw)
legend('Spectrum', 'Threshold (CFAR)', 'Detected (CFAR)', 'True', 'Fontsize', fsz)
xlabel('Range (m)', 'Fontsize', fsz)
ylabel('Amplitude (dB)', 'Fontsize', fsz)
% title('FFT result of IF')
xlim([0, range_max / 5])
% (c * 2 * pi) / (4 * pi * Ts * Slope_fre)



omegax_tau = omegalist_tau;
range_tau = (c * omegax_tau) / (4 * pi * Ts * Slope_fre)
Khat_tau = length(omegax_tau);

figure;
plot(range_idx, 10 * log10(tau_set) * ones(Nx, 1), 'r', 'Linewidth', lw, 'Markersize', msz);
hold on;
stem(range_tau, 20 * log10(abs(gainlist_tau)), 'bo', 'Linewidth', lw, 'Markersize', msz);
stem(range_true, amplitude_idx * ones(K_true, 1), ':.m', 'Linewidth', lw)
legend('Threshold (NOMP)', 'Detected (NOMP)', 'True', 'Fontsize', fsz)
xlim([0, range_max / 5])
xlabel('Range (m)', 'Fontsize', fsz)
ylabel('Amplitude (dB)', 'Fontsize', fsz)




omegax_hat = omega_list;
range_hat = (c * omegax_hat) / (4 * pi * Ts * Slope_fre)
% subplot(2, 1, 2)
figure;
plot(range_hat, 10 * log10(abs(Threshold_collect)), 'rx', 'Linewidth', lw, 'Markersize', msz);
hold on;
stem(range_hat, 20 * log10(abs(gain_list)), 'bo', 'Linewidth', lw, 'Markersize', msz);
stem(range_true, amplitude_idx * ones(K_true, 1), ':.m', 'Linewidth', lw)
legend('Threshold (NOMP-CFAR)', 'Detected (NOMP-CFAR)', 'True', 'Fontsize', fsz)
xlabel('Range (m)', 'Fontsize', fsz)
ylabel('Amplitude (dB)', 'Fontsize', fsz)
% title('NOMP-CFAR result of IF')
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




