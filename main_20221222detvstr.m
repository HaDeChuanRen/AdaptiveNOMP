% last update: 2022.7.23
% simulation of detection rate changing versus the number of observations
% in compressive scene


clc; clear; close all;

rng(5);
MC = 3000;

% Define Scenario
Nx = 256; % Length of Sinusoid
K = 16;
sigma_n = 1;              % noise variance sigma^2, instead of sigma
S_snap = 1;
SNRvec_all = 12 : 1 : 20;
length_SNR = length(SNRvec_all);

M = Nx;
Smat_com = eye(M);

% algorithm parameters
K_max = 20;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;
N_r = 60;

tau_set = sigma_n * chi2inv((1 - P_oe) ^ (1 / Nx), 2 * S_snap) / 2;
alpha_set = alpha_PoebyS(P_oe, Nx, N_r, S_snap);

% statistical variable initialize
Falsemat_tau = zeros(MC, length_SNR);
Falsemat_CA = zeros(MC, length_SNR);
Falsemat_tr = zeros(MC, length_SNR);

Detectmat_tau = zeros(MC, length_SNR);
Detectmat_CA = zeros(MC, length_SNR);
Detectmat_tr = zeros(MC, length_SNR);

% Mento Carlo method
tic;




for sp_idx = 1 : length_SNR
    SNR = SNRvec_all(sp_idx);

    for mc_idx = 1 : MC
        handle_waitbar = waitbar(((sp_idx - 1) * MC + mc_idx) / (MC * length_SNR));
        [y, omega_true, gain_true] = create_yvector(K, S_snap, SNR, sigma_n, Smat_com);

        [omegavec_tau, gainvec_tau, ~] = MNOMP(y, Smat_com, tau_set);
        resultstruct_tau = False_Detection(omega_true, gain_true,...
        omegavec_tau, gainvec_tau, Nx);
        Falsemat_tau(mc_idx, sp_idx) = resultstruct_tau.False_Eve;
        Detectmat_tau(mc_idx, sp_idx) = resultstruct_tau.Detect_Eve;

        [omegavec_CA, gainvec_CA, ~] = ...
        MNOMP_CFAR_alpha(y, Smat_com, alpha_set, N_r, K_max);
        resultstruct_CA = False_Detection(omega_true, gain_true,...
        omegavec_CA, gainvec_CA, Nx);
        Falsemat_CA(mc_idx, sp_idx) = resultstruct_CA.False_Eve;
        Detectmat_CA(mc_idx, sp_idx) = resultstruct_CA.Detect_Eve;

        [omegavec_tr, gainvec_tr, ~] = ...
        MNOMP_CFARtr(y, Smat_com, alpha_set, N_r, K_max);
        resultstruct_tr = False_Detection(omega_true, gain_true,...
        omegavec_tr, gainvec_tr, Nx);
        Falsemat_tr(mc_idx, sp_idx) = resultstruct_tr.False_Eve;
        Detectmat_tr(mc_idx, sp_idx) = resultstruct_tr.Detect_Eve;

    end
end

toc;
delete(handle_waitbar);
if MC > 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_trvsback.mat'];
    save(filename_now, 'Nx', 'P_oe', 'K', 'SNRvec_all', 'length_SNR',...
    'Falsemat_tau', 'Falsemat_CA', 'Detectmat_tau', 'Detectmat_CA',...
    'Falsemat_tr', 'Detectmat_tr');
end

% after care
Falserate_tau = mean(Falsemat_tau);
Detectrate_tau = mean(Detectmat_tau);

Falserate_CA = mean(Falsemat_CA);
Detectrate_CA = mean(Detectmat_CA);

Falserate_tr = mean(Falsemat_tr);
Detectrate_tr = mean(Detectmat_tr);





% plot the result
lw = 2;
fsz = 12;
msz = 8;


figure(1)
plot(SNRvec_all, P_oe * ones(1, length_SNR), '--k', 'Linewidth', lw)
hold on;
plot(SNRvec_all, Falserate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
plot(SNRvec_all, Falserate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
plot(SNRvec_all, Falserate_tr, '-m*', 'Linewidth', lw, 'Markersize', msz)
legend('$\bar{\rm P}_{\rm FA} = 0.01$', 'NOMP', ...
    'NOMP-CFAR', 'NOMP-CFAR (traditonal CFAR)', 'Interpreter', 'latex', ...
    'Fontsize', fsz)
xlabel('${\rm SNR}$ (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured $\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(2)
plot(SNRvec_all, Detectrate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
hold on;
plot(SNRvec_all, Detectrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
plot(SNRvec_all, Detectrate_tr, '-m*', 'Linewidth', lw, 'Markersize', msz)
xlabel('${\rm SNR}$ (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured $\bar{\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)


