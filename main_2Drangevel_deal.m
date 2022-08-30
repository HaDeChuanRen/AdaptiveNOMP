clc; clear; close all;


% read the ADC data
orginal_path = 'D:\XuMenghuai\FMCW mmwave range and Doppler estimation 202106\4program\Matlab\data\20220506exp\20220506exp';
exp_type = '\07multiple';
exp_serial = '\01';

filename = [orginal_path, exp_type, exp_serial, '\adc_data.bin'];
data_cube = readadc(filename);

range_true = [4.87, 2.63];
velocity_true = [0, 0];
amp_true = [75, 75];

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
Nx = 128;
My = 64;
NM_num = Nx * My;

N_start = 1;
M_start = 2000;
L_start = 2;
ymat = squeeze(data_cube(N_start : (N_start + Nx - 1), ...
M_start : (M_start + My - 1), L_start));

% algorithm parameter set
guard_n = 3;
guard_m = 3;
training_n = 5;
training_m = 5;
guard_training_size2D = [guard_n, guard_m, training_n, training_m];


K_max = 15;
P_oe2D = 1e-2;
N_r2D = (2 * training_n + 1) * (2 * training_m + 1) - ...
(2 * guard_n + 1) * (2 * guard_m + 1);

% NOMP method analysis
alpha_set2D = alpha_PoebyS(P_oe2D, NM_num, N_r2D);
tic;
[omegaList, gainList, ~, Threshold_collect] = ...
NOMP2D_CFAR(ymat, alpha_set2D, N_r2D, K_max);
toc;


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


lw = 2;
fsz = 12;
msz = 8;


bias = 10 * log10(NM_num);

figure;

stem3(range_hat, velocity_hat, target_threshold - bias, 'rx', 'Linewidth', lw, 'Markersize', msz);
hold on;
stem3(range_hat, velocity_hat, target_amp - bias, 'bo', 'Linewidth', lw, 'Markersize', msz);
stem3(range_true, velocity_true, amp_true - bias, ':.m', 'Linewidth', lw);

xlim([0 6])
ylim([-2 2])
xlabel('Range (m)', 'Fontsize', fsz);
ylabel('Velocity (m/s)', 'Fontsize', fsz)
zlabel('Amplitude (dB)', 'Fontsize', fsz)
legend('Threshold (NOMP-CFAR)', 'Detected (NOMP-CFAR)', 'True (people 2 and 3)', 'Fontsize', fsz)
% title('range Doppler estimation by Adap-CFAR-NOMP')





% CFAR and FFT results
tic;
ymat_absfft = abs(fft2(ymat) / sqrt(NM_num));
prob_ind_ext = repmat(ymat_absfft .^ 2, [3 3]);
cfar_detector2D = phased.CFARDetector2D('TrainingBandSize',[training_n, training_m], ...
'ThresholdFactor', 'Auto', 'GuardBandSize', [guard_n, guard_m], ...
'ProbabilityFalseAlarm', P_oe2D / NM_num, 'Method', 'CA', 'ThresholdOutputPort', true);
cfar_detector2D.NoisePowerOutputPort = true;

% figure;
% imagesc(20 * log10(ymat_absfft))

cut_idx = zeros(2, NM_num);
for n_idx = 1 : Nx
    for m_idx = 1 : My
        cut_idx(:, (n_idx - 1) * My + m_idx) = [Nx + n_idx; My + m_idx];
    end
end

[peak_grid, Threshold_CUT, sigma_hat] = cfar_detector2D(prob_ind_ext, cut_idx);
toc;
alpha_fft = mean(Threshold_CUT ./ sigma_hat);
peak_idx = find(peak_grid == 1);
peak_idx2D = cut_idx(:, peak_idx) - [Nx; My];
Khat_fft = length(peak_idx);
Threshold_fft = Threshold_CUT(peak_idx);
gain_fft = zeros(Khat_fft, 1);

for k_idx = 1 : Khat_fft
    gain_fft(k_idx) = ymat_absfft(peak_idx2D(1, k_idx), peak_idx2D(2, k_idx));
end

omegax_hatfft = 2 * pi * (peak_idx2D(1, :)' - 1) / Nx;
omegay_hatfft = wrapToPi(2 * pi * (peak_idx2D(2, :)' - 1) / My);
range_hatfft = (c * omegax_hatfft) / (4 * pi * Ts * Slope_fre);
velocity_hatfft = (c * omegay_hatfft) / (4 * pi * Fre_start * T_circle);

figure;
stem3(range_hatfft, velocity_hatfft, 10 * log10(Threshold_fft) - bias, 'rx', 'Linewidth', lw, 'Markersize', msz);
hold on;
stem3(range_hatfft, velocity_hatfft, 20 * log10(gain_fft) - bias, 'bo', 'Linewidth', lw, 'Markersize', msz);
stem3(range_true, velocity_true, amp_true - bias, ':.m', 'Linewidth', lw);
xlim([0 6])
ylim([-2 2])
xlabel('Range (m)', 'Fontsize', fsz);
ylabel('Velocity (m/s)', 'Fontsize', fsz)
zlabel('Amplitude (dB)', 'Fontsize', fsz)
legend('Threshold (CFAR)', 'Detected (CFAR)', 'True (people 2 and 3)', 'Fontsize', fsz)
% title('range Doppler estimation by FFT-CFAR')



