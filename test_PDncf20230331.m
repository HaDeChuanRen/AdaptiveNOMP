% update: 2023/3/28
% compare the calculate PD and measure PD


clc; clear; close all;
addpath('analysis tools\')
addpath('NOMP1D tools\')

rng(5);
MC = 300;

% Define Scenario
Nx = 256; % Length of Sinusoid
K = 1;
sigma_n = 1;              % noise variance sigma^2, instead of sigma
S_snap = 1;
SNRvec_all = 14 : 2 : 20;
length_SNR = length(SNRvec_all);


% algorithm parameters
P_fa = 0.01;
num_train = 30;
N_r = 2 * num_train;
num_guard = 3;

tau_set = sigma_n * chi2inv((1 - P_fa) ^ (1 / Nx), 2 * S_snap) / 2;
alpha_set = alpha_PoebyS(P_fa, Nx, N_r, S_snap);
% alpha_set = 10;

% algorithm parameters
K_max = 3;
P_fa = 0.01;
gamma_oversamping = 4;



omega_true = 2 * pi / Nx * (Nx / 2 - 0.2);
% omega_true = 2 * pi * rand();
cutvec_idx = round(Nx * omega_true / (2 * pi)) + 1;

omega_delta = (cutvec_idx - 1) * 2 * pi / Nx - omega_true + eps;
beta_omega = sin(Nx * omega_delta / 2) / (Nx * sin(omega_delta / 2));
y_full = exp(1j * (0 : (Nx - 1)).' * omega_true.') ;
% plot(abs(fft(y_full)))
train_lead = (cutvec_idx - num_guard - num_train) + (0 : num_train - 1)';
train_lead(train_lead <= 0) = train_lead(train_lead <= 0) + Nx;
train_lag = (cutvec_idx + num_guard) + (1 : num_train)';
train_lag(train_lag > Nx) = train_lag(train_lag > Nx) - Nx;
train_all = [train_lead; train_lag];

Detectmat_FFT = zeros(MC, length_SNR);
Detectmat_alpha = zeros(MC, length_SNR);
PDvec_cal = zeros(1, length_SNR);
PDvec_cal_ncf = zeros(1, length_SNR);

tic;

for sp_idx = 1 : length_SNR
    hwaitbar = waitbar((sp_idx - 1) / length_SNR);
    SNR = SNRvec_all(sp_idx);
    % fun_QT = @(T_var) exp(log(marcumq(beta_omega * sqrt(2 * 10 ^ (SNR / 10)), sqrt(2 * T_var))) + (N_r - 1) * log(T_var) - N_r * T_var / alpha_set - ...
    %     sum(log(1 : (N_r - 1))) + N_r * log(N_r / alpha_set));
    % PDvec_cal(sp_idx) = integral(fun_QT, 0, Inf);
    delta_y = 2 * (10 .^ (SNR / 10)) * (beta_omega ^ 2);
    PDvec_cal(sp_idx) = ncfcdf(alpha_set, 2, 2 * N_r, delta_y, 'upper');
    gainvec_true = sqrt(sigma_n) * (10 .^ (SNR / 20)) * exp(1j * 2 * pi * ...
        rand(K, 1)) / sqrt(Nx);
    ymat_all = gainvec_true .* y_full + sqrt(sigma_n / 2) * (randn(Nx, MC)...
        + 1j * randn(Nx, MC));

    ymat_absfft = abs(fft(ymat_all) / sqrt(Nx)) .^ 2;
    peakvec = ymat_absfft(cutvec_idx, :);
    % [peakvec, ~] = max(ymat_absfft);


    train_set = ymat_absfft(train_all, :);
    simga_vec = mean(train_set);
    Detectmat_FFT(:, sp_idx) = (peakvec > alpha_set * simga_vec)';


    for mc_idx = 1 : MC
        hwaitbar = waitbar(((sp_idx - 1) * MC + mc_idx) / (MC * length_SNR));
        gain_true = bsxfun(@times, sqrt(sigma_n) * (10 .^ (SNR / 20)), ...
            exp(1j * 2 * pi * rand(K, S_snap)));  % K* S_snap
        % y_full = exp(1j * (0 : (Nx - 1)).' * omega_true.') / sqrt(Nx) * ...
        %     gain_true;
        % y_full + sqrt(sigma_n / 2) * (randn(Nx, S_snap) + ...
        %     1j * randn(Nx, S_snap));

        y = ymat_all(:, mc_idx);
        [omegavec_alpha, gainvec_alpha, ~] = ...
            MNOMP_CFAR_alpha(y, eye(Nx), alpha_set, N_r, K_max);
        resultstruct_alpha = False_Detection(omega_true, gain_true, ...
            omegavec_alpha, gainvec_alpha, Nx);
        Detectmat_alpha(mc_idx, sp_idx) = resultstruct_alpha.Detect_Eve;
    end

end

time_MC = toc
delete(hwaitbar);

if time_MC > 600
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDcandPDm.mat'];
    save(filename_now, 'Nx', 'P_fa', 'K', 'SNRvec_all', 'length_SNR', 'MC', ...
    'PDvec_cal', 'Detectmat_FFT', 'Detectmat_alpha');
end

PDvec_alpha = squeeze(mean(Detectmat_alpha, 1));
PDvec_FFT = squeeze(mean(Detectmat_FFT, 1));

lw = 2;
fsz = 12;
msz = 8;

figure;
plot(SNRvec_all, PDvec_cal, '--k+', 'Linewidth', lw);
hold on;
% plot(SNRvec_all, PDvec_cal_ncf, '--m*', 'Linewidth', lw);
plot(SNRvec_all, PDvec_FFT, '-bo', 'Linewidth', lw, 'Markersize', msz);
plot(SNRvec_all, PDvec_alpha, '-r+', 'Linewidth', lw, 'Markersize', msz);
legend('calculated (ncf)', 'measured (FFT)', 'measured (NOMP-CFAR)', ...
    'Interpreter', 'latex', 'Fontsize', fsz);
ylabel('$\bar{\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('${\rm SNR}$ (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)



