% last update: 2022.8.29
% simulation of detection rate and false alarm rate changing
% versus the variance of noise intensity



clc; clear; close all;

rng(5);
MC = 3000;

% Define Scenario
Nx = 256; % Length of Sinusoid
K = 16;
sigma_n = 1;
S_snap = 1;
SNR = 28;

M = Nx;
Smat_com = eye(M);

% algorithm parameters
K_max = 32;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;
N_r = 60;

tau_set = sigma_n * chi2inv((1 - P_oe) ^ (1 / Nx), 2 * S_snap) / 2;
alpha_set = alpha_PoebyS(P_oe, Nx, N_r, S_snap);

u_vecall = [0; 1; 2; 4; 6; 8; 10];
length_variable = length(u_vecall);

% statistical variable initialize
Overestmat_tau = zeros(MC, length_variable);
Overestmat_CA = zeros(MC, length_variable);

Detectmat_tau = zeros(MC, length_variable);
Detectmat_CA = zeros(MC, length_variable);

% Mento Carlo method
tic;

for u_idx = 1 : length_variable
    u_set = u_vecall(u_idx);

    for mc_idx = 1 : MC
        handle_waitbar = waitbar(((u_idx - 1) * MC + mc_idx) / (MC * length_variable));
        QdB_set = (2 * rand() - 1) * u_set;
        sigma_set = 10 ^ (QdB_set / 10) * sigma_n;
        [y, omega_true, gain_true] = create_yvector(K, S_snap, SNR, sigma_set, Smat_com);

        [omegavec_tau, gainvec_tau, ~] = MNOMP(y, Smat_com, tau_set);
        resultstruct_tau = analysis_result(omega_true, gain_true,...
        omegavec_tau, gainvec_tau, Nx, gamma_oversamping);
        Overestmat_tau(mc_idx, u_idx) = resultstruct_tau.Overest_Eve;
        Detectmat_tau(mc_idx, u_idx) = resultstruct_tau.Detect_Eve;

        [omegavec_CA, gainvec_CA, ~] = ...
        MNOMP_CFAR_alpha(y, Smat_com, alpha_set, N_r, K_max);
        resultstruct_CA = analysis_result(omega_true, gain_true,...
        omegavec_CA, gainvec_CA, Nx, gamma_oversamping);
        Overestmat_CA(mc_idx, u_idx) = resultstruct_CA.Overest_Eve;
        Detectmat_CA(mc_idx, u_idx) = resultstruct_CA.Detect_Eve;
    end
end


toc;
delete(handle_waitbar);

% after care
Overestrate_tau = mean(Overestmat_tau);
Detectrate_tau = mean(Detectmat_tau);

Overestrate_CA = mean(Overestmat_CA);
Detectrate_CA = mean(Detectmat_CA);

if MC > 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDvsSNR.mat'];
    save(filename_now, 'Nx', 'P_oe', 'K', 'u_vecall', 'length_SNR',...
    'Overestmat_tau', 'Overestmat_CA', 'Detectmat_tau', 'Detectmat_CA');
end

% plot the result
lw = 1.6;
fsz = 12;
msz = 10;


figure(1)
plot(u_vecall, P_oe * ones(1, length_variable), '--k', 'Linewidth', lw)
hold on;
plot(u_vecall, Overestrate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
plot(u_vecall, Overestrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
legend('${\rm P}_{\rm OE} = 0.01$', 'NOMP', ...
    'NOMP-CA', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('${\rm SNR}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(2)
plot(u_vecall, Detectrate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
hold on;
plot(u_vecall, Detectrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
xlabel('${\rm SNR}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured ${\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)
