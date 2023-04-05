clc; clear; close all;


% read the ADC data

orginal_path = 'C:\study\202206MNOMP_CFAR\data\20220506exp';

exp_type = '\02people2';
exp_serial = '\01';
N_start = 1;
M_start = 17; %17 1850
% L_start = 1;

range_true = [2.95 + 0.14; 4.78 + 0.14];
theta_true = [- 19.8; 0];
K_ture = length(range_true);
amp_true = 60 * ones(K_ture, 1);

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


ymat = squeeze(data_cube(N_start : (N_start + Nx - 1), M_start, :));

% algorithm parameter set
guard_n = 4;
guard_l = 0;
training_n = 12;
training_l = 1;
guard_training_size2D = [guard_n, guard_l, training_n, training_l];


K_max = 32;
P_oe2D = 4e-4;
N_r2D = (2 * training_n + 1) * (2 * training_l + 1) - ...
(2 * guard_n + 1) * (2 * guard_l + 1);
% N_r2D = 40;
alpha_set2D = alpha_PoebyS(P_oe2D, NL_num, N_r2D);
guard_band = [guard_n, guard_l];
sigma_set = 10 ^ (24 / 10);
tau_set = sigma_set * chi2inv((1 - P_oe2D) ^ (1 / NL_num), 2) / 2;

% NOMP method analysis
tic;
[omega_NOMPCFAR, gain_NOMPCFAR, ~, Threshold_NOMPCFAR] = ...
NOMP2D_fast(ymat, alpha_set2D, N_r2D, K_max, guard_band);
time_NOMPCFAR = toc;

tic;
[omega_tau, gain_tau, ~] = ...
NOMP2D(ymat, tau_set, K_max);
time_tau = toc;

omegax_NOMPCFAR = omega_NOMPCFAR(:, 1);
omegaz_NOMPCFAR = wrapToPi(omega_NOMPCFAR(:, 2));
range_NOMPCFAR = (c * omegax_NOMPCFAR) / (4 * pi * Ts * Slope_fre);
theta_NOMPCFAR = asin((c * omegaz_NOMPCFAR) / (2 * pi * Fre_start * Rx_interval));
theta_NOMPCFAR_deg = theta_NOMPCFAR * 180 / pi;



lw = 2;
fsz = 12;
msz = 8;
range_max = (c * 2 * pi) / (4 * pi * Ts * Slope_fre);
% bias = 10 * log10(NL_num);

figure;
plot3(range_NOMPCFAR, theta_NOMPCFAR_deg, 10 * log10(abs(Threshold_NOMPCFAR)), 'rx', 'Linewidth', lw, 'Markersize', msz);
hold on;
stem3(range_NOMPCFAR, theta_NOMPCFAR_deg, 20 * log10(abs(gain_NOMPCFAR)), 'bo', 'Linewidth', lw, 'Markersize', msz);
stem3(range_true, theta_true, amp_true, ':.m', 'Linewidth', lw); 
grid on;
xlim([0 range_max / 2])
ylim([-40 40])
xlabel('Range (m)', 'Fontsize', fsz);
ylabel('Azimuth ($\circ$)', 'Interpreter','latex', 'Fontsize', fsz);
zlabel('Amplitude (dB)', 'Fontsize', fsz)
legend('Threshold (NOMP-CFAR)', 'Amplitude (NOMP-CFAR)', 'True', 'Fontsize', fsz)
% title('location estimation by Adap-CFAR-NOMP')

omegax_tau = omega_tau(:, 1);
omegaz_tau = wrapToPi(omega_tau(:, 2));
range_tau = (c * omegax_tau) / (4 * pi * Ts * Slope_fre);
theta_tau = asin((c * omegaz_tau) / (2 * pi * Fre_start * Rx_interval));
theta_tau_deg = theta_tau * 180 / pi;
Khat_tau = length(omegax_tau);


range_idx = linspace(0, range_max, Nx);
theta_idx = linspace(-90, 90, Lz);
[range_newidx, theta_newidx] = meshgrid(range_idx, theta_idx);

figure;
surf(range_newidx, theta_newidx, (10 * log10(tau_set)) * ones(Lz, Nx),...
'Linestyle', 'none');
hold on;
stem3(range_tau, theta_tau_deg, 20 * log10(abs(gain_tau)), ...
'bo', 'Linewidth', lw, 'Markersize', msz);
stem3(range_true, theta_true, amp_true, ':.m', 'Linewidth', lw); 'True',
grid on;
xlim([0 range_max / 2])
ylim([-40 40])
xlabel('Range (m)', 'Fontsize', fsz);
ylabel('Azimuth ($\circ$)', 'Interpreter','latex', 'Fontsize', fsz);
zlabel('Amplitude (dB)', 'Fontsize', fsz)
legend('Threshold (NOMP)', 'Amplitude (NOMP)', 'Fontsize', fsz)


% CFAR and FFT results
% tic;
% ymat_absfft = abs(fft2(ymat) / sqrt(NL_num));
% prob_ind_ext = repmat(ymat_absfft .^ 2, [3 3]);
% cfar_detector2D = phased.CFARDetector2D('TrainingBandSize',[training_n, training_l], ...
% 'ThresholdFactor', 'Auto', 'GuardBandSize', [guard_n, guard_l], ...
% 'ProbabilityFalseAlarm', P_oe2D / NL_num, 'Method', 'CA', 'ThresholdOutputPort', true);
% cfar_detector2D.NoisePowerOutputPort = true;
% 
% % figure;
% % imagesc(20 * log10(ymat_absfft))
% 
% cut_idx = zeros(2, NL_num);
% for n_idx = 1 : Nx
%     for l_idx = 1 : Lz
%         cut_idx(:, (n_idx - 1) * Lz + l_idx) = [Nx + n_idx; Lz + l_idx];
%     end
% end
% 
% [peak_grid, Threshold_CUT, sigma_hat] = cfar_detector2D(prob_ind_ext, cut_idx);
% toc;
% alpha_fft = mean(Threshold_CUT ./ sigma_hat);
% peak_idx = find(peak_grid == 1);
% peak_idx2D = cut_idx(:, peak_idx) - [Nx; Lz];
% Khat_fft = length(peak_idx);
% Threshold_fft = Threshold_CUT(peak_idx);
% gain_fft = zeros(Khat_fft, 1);
% 
% for k_idx = 1 : Khat_fft
%     gain_fft(k_idx) = ymat_absfft(peak_idx2D(1, k_idx), peak_idx2D(2, k_idx));
% end
% 
% omegax_hatfft = 2 * pi * (peak_idx2D(1, :)' - 1) / Nx;
% omegaz_hatfft = wrapToPi(2 * pi * (peak_idx2D(2, :)' - 1) / Lz);
% range_hatfft = (c * omegax_hatfft) / (4 * pi * Ts * Slope_fre);
% theta_hatfft = asin((c * omegaz_hatfft) / (2 * pi * Fre_start * Rx_interval));
% 
% xloc_hatfft = range_hatfft .* sin(theta_hatfft);
% yloc_hatfft = range_hatfft .* cos(theta_hatfft);
% 
% figure;
% stem3(xloc_hatfft, yloc_hatfft, 20 * log10(abs(gain_fft)));
% hold on;
% plot3(xloc_hatfft, yloc_hatfft, 10 * log10(abs(Threshold_fft)), 'r*');
% xlim([-5 5])
% ylim([0 6])
% xlabel('xaxis(m)');
% ylabel('yaxis(m)');
% zlabel('amplitude(dB)')
% legend('targets amplitude', 'thereshold')
% title('location estimation by FFT-CFAR')


