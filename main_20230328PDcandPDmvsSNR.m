% update: 2023/3/28
% compare the calculate PD and measure PD


clc; clear; close all;
addpath('analysis tools\')
addpath('NOMP1D tools\')

rng(5);
MC = 300;

% Define Scenario
Nx = 1024; % Length of Sinusoid
K = 1;
sigma_n = 1;              % noise variance sigma^2, instead of sigma
S_snap = 1;
SNRvec_all = 6 : 2 : 20;
length_SNR = length(SNRvec_all);

My = Nx;
Smat_com = eye(My);

% algorithm parameters
K_max = 3;
P_fa = 0.01;
gamma_oversamping = 4;
num_train = 100;
N_r = 2 * num_train;
num_guard = 4;

tau_set = sigma_n * chi2inv((1 - P_fa) ^ (1 / Nx), 2 * S_snap) / 2;
alpha_set = alpha_PoebyS(P_fa, Nx, N_r, S_snap);


delta_all = [0; 0.25; 0.49];

length_delta = length(delta_all);
% length_delta = 1;
% omega_true = 2 * pi / Nx * (Nx / 2);
% omega_delta = (cutvec_idx - 1) * 2 * pi / Nx - omega_true + eps;

% cutvec_idx = round(Nx * omega_true / (2 * pi)) + 1;
% cutvec_tar = cutvec_idx * ones(1, MC);
cutvec_idx = Nx / 2 - 5;
trainmat_lead = (cutvec_idx - num_guard - num_train) + (0 : num_train - 1)';
trainmat_lead(trainmat_lead <= 0) = trainmat_lead(trainmat_lead <= 0) + Nx;
trainmat_lag = (cutvec_idx + num_guard) + (1 : num_train)';
trainmat_lag(trainmat_lag > Nx) = trainmat_lag(trainmat_lag > Nx) - Nx;
train_mat = [trainmat_lead; trainmat_lag];

Detectmat_FFT = zeros(MC, length_SNR, length_delta);
Detectmat_alpha = zeros(MC, length_SNR, length_delta);
PDvec_cal = zeros(length_SNR, length_delta);


tic;
for d_idx = 1 : length_delta

    delta_d = delta_all(d_idx);
    omega_delta = delta_d * 2 * pi / Nx + eps;
    omega_true = (cutvec_idx + delta_d - 1) * 2 * pi / Nx;
    beta_omega = sin(Nx * omega_delta / 2) / (Nx * sin(omega_delta / 2));
    y_full = exp(1j * (0 : (Nx - 1)).' * omega_true.') / sqrt(Nx);
    ymat_full = repmat(y_full, [1, MC]);

    for sp_idx = 1 : length_SNR
        hwaitbar = waitbar(((d_idx - 1) * length_SNR + sp_idx) / ...
            (length_delta * length_SNR));

        SNR = SNRvec_all(sp_idx);
        delta_y = 2 * (10 .^ (SNR / 10)) * (beta_omega ^ 2);
        PDvec_cal(sp_idx, d_idx) = ncfcdf(alpha_set, 2, 2 * N_r, delta_y, ...
            'upper');
        gainvec_true = sqrt(sigma_n) * (10 .^ (SNR / 20)) * exp(1j * 2 * ...
            pi * rand(K, MC));
        ymat_all = gainvec_true .* ymat_full + sqrt(sigma_n / 2) * ...
            (randn(Nx, MC) + 1j * randn(Nx, MC));

        % debug---------
        % if (sp_idx ~= 4) || (d_idx ~= 3)
        %     continue;
        % end
        % mc_idx = 1193;
        % peakvec(mc_idx)
        % sigmavec(mc_idx)
        % ----------------


        ymat_absfft = abs(fft(ymat_all) / sqrt(Nx)) .^ 2;
        % plot(ymat_absfft(:, 1))
        peakvec = ymat_absfft(cutvec_idx, :);

        train_set = ymat_absfft(train_mat, :);
        sigmavec = mean(train_set);
        Detectmat_FFT(:, sp_idx, d_idx) = (peakvec > alpha_set * sigmavec)';

        for mc_idx = 1 : MC
            y = ymat_all(:, mc_idx);
            gain_true = gainvec_true(:, mc_idx);
            [omegavec_alpha, gainvec_alpha, ~] = ...
                MNOMP_CFAR_alpha(y, Smat_com, alpha_set, N_r, K_max);
            if ~isempty(omegavec_alpha) && min(abs(omegavec_alpha - ...
                omega_true)) < (pi / Nx)
                Detectmat_alpha(mc_idx, sp_idx, d_idx) = 1;
            end
        end
    end
end

time_MC = toc
delete(hwaitbar);

% debug----------
% dismat_NF = squeeze(Detectmat_alpha(:, :, 3) - Detectmat_FFT(:, :, 3));
% disidx_NF = find(dismat_NF == 1);
% [disrow_NF, discol_NF] = ind2sub([MC, length_SNR], disidx_NF);
% -----------------


if time_MC > 600
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDcandPDm.mat'];
    save(filename_now, 'Nx', 'P_fa', 'K', 'SNRvec_all', 'length_SNR', 'MC',...
    'PDvec_cal', 'Detectmat_FFT', 'Detectmat_alpha');
end

PDvec_alpha = squeeze(mean(Detectmat_alpha, 1));
PDvec_FFT = squeeze(mean(Detectmat_FFT, 1));

lw = 2;
fsz = 10;
msz = 8;

figure;
hold on;

plot(SNRvec_all, PDvec_cal(:, 1), '--k+', 'Linewidth', lw, ...
    'Markersize', msz);
plot(SNRvec_all, PDvec_FFT(:, 1), '-b*', 'Linewidth', lw, ...
    'Markersize', msz);
plot(SNRvec_all, PDvec_alpha(:, 1), '-ro', 'Linewidth', lw, 'Markersize', msz);

plot(SNRvec_all, PDvec_cal(:, 2), '--kx', 'Linewidth', lw, ...
    'Markersize', msz);
plot(SNRvec_all, PDvec_FFT(:, 2), '-b+', 'Linewidth', lw, ...
    'Markersize', msz);
plot(SNRvec_all, PDvec_alpha(:, 2), '-r^', 'Linewidth', lw, 'Markersize', msz);

plot(SNRvec_all, PDvec_cal(:, 3), '--kd', 'Linewidth', lw, ...
    'Markersize', msz);
plot(SNRvec_all, PDvec_FFT(:, 3), '-bo', 'Linewidth', lw, ...
    'Markersize', msz);
plot(SNRvec_all, PDvec_alpha(:, 3), '-rv', 'Linewidth', lw, 'Markersize', msz);

legend('eq.(41) ($\delta = 0$)', 'FFT ($\delta = 0$)', ...
'NOMP-CFAR ($\delta = 0$)', 'eq.(41) ($\delta = 0.25$)', ...
'FFT ($\delta = 0.25$)', 'NOMP-CFAR ($\delta = 0.25$)', ...
'eq.(41) ($\delta = 0.49$)', 'FFT ($\delta = 0.49$)', ...
'NOMP-CFAR ($\delta = 0.49$)', 'Interpreter', 'latex', 'Fontsize', fsz);
ylabel('$\bar{\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('${\rm SNR}$ (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)



