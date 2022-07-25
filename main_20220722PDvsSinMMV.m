% last update: 2022.7.22
% simulation of detection rate chenging versus the number of snapshots


clc; clear; close all;

rng(5);

MC = 50;

% Define Scenario
Nx = 256; % Length of Sinusoid
K = 16;
sigma_n = 1;              % noise variance sigma^2, instead of sigma


% Svec_all = [1, 3, 5, 8, 10, 20, 30, 40, 50];
Svec_all = 1 : 8;
length_S = length(Svec_all);
SNR = 12;

M = Nx;
Smat_com = eye(M);

% algorithm parameters
K_max = 20;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;
N_r = 60;

% statistical variable initialize
Overestmat_tau = zeros(MC, length_S);
Overestmat_CA = zeros(MC, length_S);

Detectmat_tau = zeros(MC, length_S);
Detectmat_CA = zeros(MC, length_S);


% Mento Carlo method
tic;

for sp_idx = 1 : length_S
    S_snap = Svec_all(sp_idx);
    alpha_set = alpha_PoebyS(P_oe, Nx, N_r, S_snap);
    tau_set = sigma_n * chi2inv((1 - P_oe) ^ (1 / Nx), 2 * S_snap) / 2;

    for mc_idx = 1 : MC
        handle_waitbar = waitbar(((sp_idx - 1) * MC + mc_idx) / (MC * length_S));
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
plot(Svec_all(1 : 6), P_oe * ones(1, length_S - 3), '--k', 'Linewidth', lw)
hold on;
plot(Svec_all(1 : 6), Overestrate_tau(1 : 6), '-ro', 'Linewidth', lw, 'Markersize', msz)
plot(Svec_all(1 : 6), Overestrate_CA(1 : 6), '-b+', 'Linewidth', lw, 'Markersize', msz)
legend('${\rm P}_{\rm OE} = 0.01$', 'NOMP', ...
    'NOMP-CA', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('Number of snapshot $S$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(2)
plot(Svec_all(1 : 6), Detectrate_tau(1 : 6), '-ro', 'Linewidth', lw, 'Markersize', msz)
hold on;
plot(Svec_all(1 : 6), Detectrate_CA(1 : 6), '-b+', 'Linewidth', lw, 'Markersize', msz)
xlabel('Number of snapshots $S$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured ${\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)

if MC > 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDvsSNR.mat'];
    save(filename_now, 'Nx', 'P_oe', 'K', 'Svec_all', 'length_S',...
    'Overestmat_tau', 'Overestmat_CA', 'Detectmat_tau', 'Detectmat_CA');
end


