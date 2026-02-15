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
SNR = 22;

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

u_vecall = [0; 2; 4; 6; 8; 10];
length_variable = length(u_vecall);

% statistical variable initialize
Falsemat_tau = zeros(MC, length_variable);
Detectmat_tau = zeros(MC, length_variable);

Falsemat_CA = zeros(MC, length_variable);
Detectmat_CA = zeros(MC, length_variable);

Falsemat_VALSE = zeros(MC, length_variable);
Detectmat_VALSE = zeros(MC, length_variable);

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
        resultstruct_tau = False_Detection(omega_true, gain_true,...
        omegavec_tau, gainvec_tau, Nx);
        Falsemat_tau(mc_idx, u_idx) = resultstruct_tau.False_Eve;
        Detectmat_tau(mc_idx, u_idx) = resultstruct_tau.Detect_Eve;

        [omegavec_CA, gainvec_CA, ~] = ...
        MNOMP_CFAR_alpha(y, Smat_com, alpha_set, N_r, K_max);
        resultstruct_CA = False_Detection(omega_true, gain_true,...
        omegavec_CA, gainvec_CA, Nx);
        Falsemat_CA(mc_idx, u_idx) = resultstruct_CA.False_Eve;
        Detectmat_CA(mc_idx, u_idx) = resultstruct_CA.Detect_Eve;

        y_full = exp(1j * (0:(Nx-1)).' * omega_true.') / sqrt(Nx) * gain_true;
        out_VALSE = MVALSE_best(y, (0 : (Nx - 1))', 2, y_full);
        omegaList_VALSE = wrapTo2Pi(out_VALSE.freqs);
        gainList_VALSE = out_VALSE.amps;
        resultstruct_VALSE = False_Detection(omega_true, gain_true,...
        omegaList_VALSE, gainList_VALSE, Nx);
        Falsemat_VALSE(mc_idx, u_idx) = resultstruct_VALSE.False_Eve;
        Detectmat_VALSE(mc_idx, u_idx) = resultstruct_VALSE.Detect_Eve;
    end
end


toc;
delete(handle_waitbar);

% after care
Falserate_tau = mean(Falsemat_tau);
Detectrate_tau = mean(Detectmat_tau);

Falserate_CA = mean(Falsemat_CA);
Detectrate_CA = mean(Detectmat_CA);

Falserate_VALSE = mean(Falsemat_VALSE);
Detectrate_VALSE = mean(Detectmat_VALSE);

if MC > 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDvsSNR.mat'];
    save(filename_now, 'Nx', 'P_oe', 'K', 'u_vecall', 'length_variable',...
    'Falsemat_tau', 'Falsemat_CA', 'Detectmat_tau', 'Detectmat_CA', ...
    'Falsemat_VALSE', 'Detectmat_VALSE');
end

% plot the result
lw = 2;
fsz = 12;
msz = 8;


figure(1)
semilogy(u_vecall, P_oe * ones(1, length_variable) / Nx, '--k', 'Linewidth', lw)
hold on;
semilogy(u_vecall, Falserate_VALSE / Nx, '-md', 'Linewidth', lw, 'Markersize', msz)
semilogy(u_vecall, Falserate_tau / Nx, '-ro', 'Linewidth', lw, 'Markersize', msz)
semilogy(u_vecall, Falserate_CA / Nx, '-b+', 'Linewidth', lw, 'Markersize', msz)
legend('$\bar{\rm P}_{\rm FA} = 0.01 / N$', 'VALSE', 'NOMP', ...
    'NOMP-CFAR', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('Strength of Noise fluctuation u (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured $\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)


figure(2)
hold on;
plot(u_vecall, Detectrate_VALSE, '-md', 'Linewidth', lw, 'Markersize', msz)
plot(u_vecall, Detectrate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
plot(u_vecall, Detectrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
xlabel('Strength of Noise Fluctuation u (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured ${\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(3)
plot(u_vecall, P_oe * ones(1, length_variable) / Nx, '--k', 'Linewidth', lw)
hold on;
plot(u_vecall, Falserate_VALSE / Nx, '-md', 'Linewidth', lw, 'Markersize', msz)
plot(u_vecall, Falserate_tau / Nx, '-ro', 'Linewidth', lw, 'Markersize', msz)
plot(u_vecall, Falserate_CA / Nx, '-b+', 'Linewidth', lw, 'Markersize', msz)
legend('$\bar{\rm P}_{\rm FA} = 0.01 / N$', 'VALSE', 'NOMP', ...
    'NOMP-CFAR', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('Strength of Noise fluctuation u (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured $\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylim([0, 1e-3])


