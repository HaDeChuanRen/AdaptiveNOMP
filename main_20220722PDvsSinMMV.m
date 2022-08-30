% last update: 2022.7.22
% simulation of detection rate chenging versus the number of snapshots


clc; clear; close all;

rng(5);

MC = 3000;

% Define Scenario
Nx = 256; % Length of Sinusoid
K = 16;
sigma_n = 1;              % noise variance sigma^2, instead of sigma


% Svec_all = [1, 3, 5, 8, 10, 20, 30, 40, 50];
Svec_all = 1 : 8;
length_S = length(Svec_all);
SNR = 10;

M = Nx;
Smat_com = eye(M);

% algorithm parameters
K_max = 20;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;
N_r = 60;

% statistical variable initialize
Falsemat_tau = zeros(MC, length_S);
Falsemat_CA = zeros(MC, length_S);

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
        resultstruct_tau = False_Detection(omega_true, gain_true,...
        omegavec_tau, gainvec_tau, Nx);
        Falsemat_tau(mc_idx, sp_idx) = resultstruct_tau.False_Eve;
        Detectmat_tau(mc_idx, sp_idx) = resultstruct_tau.Detect_Eve;

        [omegavec_CA, gainvec_CA, ~] = ...
        MNOMP_forward_alpha(y, Smat_com, alpha_set, N_r, K_max);
        resultstruct_CA = False_Detection(omega_true, gain_true,...
        omegavec_CA, gainvec_CA, Nx);
        Falsemat_CA(mc_idx, sp_idx) = resultstruct_CA.False_Eve;
        Detectmat_CA(mc_idx, sp_idx) = resultstruct_CA.Detect_Eve;

    end
end

toc;
delete(handle_waitbar);

% after care
Falserate_tau = mean(Falsemat_tau);
Detectrate_tau = mean(Detectmat_tau);

Falserate_CA = mean(Falsemat_CA);
Detectrate_CA = mean(Detectmat_CA);

% plot the result
lw = 2;
fsz = 12;
msz = 8;


figure(1)
plot(Svec_all, P_oe * ones(1, length_S), '--k', 'Linewidth', lw)
hold on;
plot(Svec_all, Falserate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
plot(Svec_all, Falserate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
legend('$\bar{\rm P}_{\rm FA} = 0.01$', 'NOMP', ...
    'NOMP-CFAR', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('Number of snapshot $S$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured $\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(2)
plot(Svec_all, Detectrate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
hold on;
plot(Svec_all, Detectrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
xlabel('Number of snapshots $S$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured ${\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)

if MC > 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDinMMV.mat'];
    save(filename_now, 'Nx', 'P_oe', 'K', 'Svec_all', 'length_S',...
    'Falsemat_tau', 'Falsemat_CA', 'Detectmat_tau', 'Detectmat_CA');
end


