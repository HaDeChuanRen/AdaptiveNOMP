% last update: 2023/4/3

clc; clear; close all;

addpath('analysis tools\')
addpath('NOMP1D tools\')

rng(4);
MC = 3000;

% Define Scenario
Nx = 256; % Length of Sinusoid
K = 16;
sigma_n = 1;              % noise variance sigma^2, instead of sigma
S_snap = 1;
SNRvec_all = 6 : 2 : 20;
SNR_ord = 20;
length_SNR = length(SNRvec_all);
omega_min = 2 * pi / Nx;

My = Nx;
Smat_com = eye(My);

cutvec_idx = Nx / 2 - 5;
delta_d = 0.2;
omega_delta = delta_d * 2 * pi / Nx + eps;
omega_det = (cutvec_idx + delta_d - 1) * 2 * pi / Nx;
beta_omega = sin(Nx * omega_delta / 2) / (Nx * sin(omega_delta / 2));
yvec_det = exp(1j * (0 : (Nx - 1)).' * omega_det.') / sqrt(Nx);

% NOMP-CFAR algorithm parameters
K_max = 32;
Pfa_low = 0.01;
Pfa_high = 0.05;
gamma_oversamping = 4;
num_train = 30;
N_r = 2 * num_train;
num_guard = 4;
alpha_low = alpha_PoebyS(Pfa_low, Nx, N_r, S_snap);
alpha_high = alpha_PoebyS(Pfa_high, Nx, N_r, S_snap);
% tau_set = sigma_n * chi2inv((1 - P_fa) ^ (1 / Nx), 2 * S_snap) / 2;

% FFT parameters
trainmat_lead = (cutvec_idx - num_guard - num_train) + (0 : num_train - 1)';
trainmat_lead(trainmat_lead <= 0) = trainmat_lead(trainmat_lead <= 0) + Nx;
trainmat_lag = (cutvec_idx + num_guard) + (1 : num_train)';
trainmat_lag(trainmat_lag > Nx) = trainmat_lag(trainmat_lag > Nx) - Nx;
train_mat = [trainmat_lead; trainmat_lag];


% initialization
Detmat_alphalow = zeros(MC, length_SNR);
Detmat_FFTlow = zeros(MC, length_SNR);
PDcal_low = zeros(1, length_SNR);

Detmat_alphahigh = zeros(MC, length_SNR);
Detmat_FFThigh = zeros(MC, length_SNR);
PDcal_high = zeros(1, length_SNR);


tic;
for sp_idx = 1 : length_SNR
    hwaitbar = waitbar((sp_idx - 1) / length_SNR);

    SNR_det = SNRvec_all(sp_idx);
    SNR = [SNR_det; SNR_ord * ones(K - 1, 1)];
    delta_y = 2 * (10 .^ (SNR_det / 10)) * (beta_omega ^ 2);
    PDcal_low(sp_idx) = ncfcdf(alpha_low, 2, 2 * N_r, delta_y, 'upper');
    PDcal_high(sp_idx) = ncfcdf(alpha_high, 2, 2 * N_r, delta_y, 'upper');

    for mc_idx = 1 : MC
        omega_true = zeros(K, 1);
        omega_true(1) = omega_det;
        for k = 2 : K
            th = pi * (2 * rand - 1);
            while min(abs((wrapToPi(th - omega_true(1 : k - 1))))) < omega_min
                th = pi * (2 * rand() - 1);
            end
            omega_true(k) = th;
        end
        omega_true = wrapTo2Pi(omega_true);
        gain_true = sqrt(sigma_n) * (10 .^ (SNR / 20)) .* exp(1j * 2 * pi * ...
            rand(K, S_snap));

        yfull_det = gain_true(1, :) * yvec_det;
        y_noise = sqrt(sigma_n / 2) * (randn(Nx, 1) + 1j * randn(Nx, 1));

        ydet = yfull_det + y_noise;
        y_full =  exp(1j * (0 : (Nx - 1)).' * omega_true.') / sqrt(Nx) * ...
            gain_true;
        y = y_full + y_noise;

        ydet_absfft = abs(fft(ydet) / sqrt(Nx)) .^ 2;
        % plot(ymat_absfft)
        peakvec = ydet_absfft(cutvec_idx, :);
        train_set = ydet_absfft(train_mat, :);
        sigmavec = mean(train_set);
        Detmat_FFTlow(mc_idx, sp_idx) = (peakvec > alpha_low * sigmavec);
        Detmat_FFThigh(mc_idx, sp_idx) = (peakvec > alpha_high * sigmavec);

        [omega_alphalow, ~, ~] = MNOMP_CFAR_alpha(y, Smat_com, alpha_low, ...
            N_r, K_max);
        if ~isempty(omega_alphalow) && min(abs(omega_alphalow - omega_det))...
            < (pi / Nx)
            Detmat_alphalow(mc_idx, sp_idx) = 1;
        end

        [omega_alphahigh, ~, ~] = MNOMP_CFAR_alpha(y, Smat_com, alpha_high, ...
            N_r, K_max);
        if ~isempty(omega_alphahigh) && min(abs(omega_alphahigh - omega_det))...
            < (pi / Nx)
            Detmat_alphahigh(mc_idx, sp_idx) = 1;
        end
    end
end

time_MC = toc;
delete(hwaitbar);

if time_MC > 600 && MC > 1000
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDcnvsSNRmul.mat'];
    save(filename_now, 'Nx', 'Pfa_low', 'Pfa_high', 'K', 'SNRvec_all', ...
        'length_SNR', 'MC', 'PDcal_low', 'Detmat_alphalow', 'PDcal_high', ...
        'Detmat_alphahigh');
end

PDmea_alphalow = mean(Detmat_alphalow, 1);
PDmea_alphahigh = mean(Detmat_alphahigh, 1);
PDmea_FFTlow = mean(Detmat_FFTlow, 1);
PDmea_FFThigh = mean(Detmat_FFThigh, 1);


lw = 2;
fsz = 12;
msz = 8;

figure;
hold on;

plot(SNRvec_all, PDcal_low, '--k+', 'Linewidth', lw, 'Markersize', msz);
plot(SNRvec_all, PDmea_FFTlow, '-bo', 'Linewidth', lw, 'Markersize', msz);
plot(SNRvec_all, PDmea_alphalow, '-r^', 'Linewidth', lw, 'Markersize', msz);

plot(SNRvec_all, PDcal_high, '--kx', 'Linewidth', lw, 'Markersize', msz);
plot(SNRvec_all, PDmea_FFThigh, '-bd', 'Linewidth', lw, 'Markersize', msz);
plot(SNRvec_all, PDmea_alphahigh, '-ro', 'Linewidth', lw, 'Markersize', msz);

legend('eq.(41) ($\bar{\rm P}_{\rm FA} = 0.01$)', ...
'FFT (oracle) ($\bar{\rm P}_{\rm FA} = 0.01$)', ...
'NOMP-CFAR ($\bar{\rm P}_{\rm FA} = 0.01$)', ...
'eq.(41) ($\bar{\rm P}_{\rm FA} = 0.05$)', ...
'FFT (oracle) ($\bar{\rm P}_{\rm FA} = 0.05$)', ...
'NOMP-CFAR ($\bar{\rm P}_{\rm FA} = 0.05$)', ...
'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('$\bar{\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('${\rm SNR}$ (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)

