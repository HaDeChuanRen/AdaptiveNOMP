% last update: 2022.7.23
% simulation of detection rate changing versus the number of observations
% in compressive scene


clc; clear; close all;

rng(5);
MC = 5000;

% Define Scenario
Nx = 256; % Length of Sinusoid
K = 8;
sigma_n = 1;              % noise variance sigma^2, instead of sigma
S_snap = 1;
SNR = 22;

MNvec_ratio_all = (1 : 8) / 8;
length_MNratio = length(MNvec_ratio_all);
% type of measurement matrix
measure_type = 'cmplx_bernoulli';
% windowing weights
% options 'all_ones' (default), 'hamming' and 'hann'
window_type = 'all_ones';

% algorithm parameters
K_max = K * 2;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;
N_r = 60;
alpha_set = alpha_PoebyS(P_oe, Nx, N_r, S_snap);
tau_set = sigma_n * chi2inv((1 - P_oe) ^ (1 / Nx), 2 * S_snap) / 2;

% statistical variable initialize
Overestmat_tau = zeros(MC, length_MNratio);
Overestmat_CA = zeros(MC, length_MNratio);

Detectmat_tau = zeros(MC, length_MNratio);
Detectmat_CA = zeros(MC, length_MNratio);

% Mento-Carlo method
tic;
for sp_idx = 1 : length_MNratio
    M = Nx * MNvec_ratio_all(sp_idx);
    Smat_com = generateMeasMat(Nx, M, measure_type, window_type);

    for mc_idx = 1 : MC
        handle_waitbar = waitbar(((sp_idx - 1) * MC + mc_idx) / (MC * length_MNratio));
        [y, omega_true, gain_true] = create_yvector(K, S_snap, SNR, sigma_n, Smat_com);

        [omegavec_tau, gainvec_tau, ~] = MNOMP(y, Smat_com, tau_set);
        resultstruct_tau = analysis_result(omega_true, gain_true,...
        omegavec_tau, gainvec_tau, Nx, gamma_oversamping);
        Overestmat_tau(mc_idx, sp_idx) = resultstruct_tau.Overest_Eve;
        Detectmat_tau(mc_idx, sp_idx) = resultstruct_tau.Detect_Eve;

        [omegavec_CA, gainvec_CA, ~] = ...
        MNOMP_CFAR_alpha(y, Smat_com, alpha_set, N_r, K_max);
        resultstruct_CA = analysis_result(omega_true, gain_true,...
        omegavec_CA, gainvec_CA, Nx, gamma_oversamping);
        Overestmat_CA(mc_idx, sp_idx) = resultstruct_CA.Overest_Eve;
        Detectmat_CA(mc_idx, sp_idx) = resultstruct_CA.Detect_Eve;

    end
end

toc;
delete(handle_waitbar);

% after care
Overestrate_tau = mean(Overestmat_tau);
Detectrate_tau = mean(Detectmat_tau);

Overestrate_CA = mean(Overestmat_CA);
Detectrate_CA = mean(Detectmat_CA);

% plot the result
lw = 1.6;
fsz = 12;
msz = 10;


figure(1)
plot(MNvec_ratio_all, P_oe * ones(1, length_MNratio), '--k', 'Linewidth', lw)
hold on;
plot(MNvec_ratio_all, Overestrate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
plot(MNvec_ratio_all, Overestrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
legend('${\rm P}_{\rm OE} = 0.01$', 'NOMP', ...
    'NOMP-CA', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('Compressive rate', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(2)
plot(MNvec_ratio_all, Detectrate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
hold on;
plot(MNvec_ratio_all, Detectrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
xlabel('Compressive rate', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured ${\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)



if MC > 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDvsSNR.mat'];
    save(filename_now, 'Nx', 'P_oe', 'K', 'MNvec_ratio_all', 'length_MNratio',...
    'Overestmat_tau', 'Overestmat_CA', 'Detectmat_tau', 'Detectmat_CA');
end

