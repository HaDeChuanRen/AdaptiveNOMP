% last update: 2022.12.12
% simulation of detection rate and false alarm rate changing
% versus the variance of noise intensity



clc; clear; close all;

rng(8);
MC = 2000;

% Define Scenario
Nx = 256; % Length of Sinusoid
K = 16;
sigma_n = 1;
S_snap = 1;
SNR = 24;
mu_sigma = 1;

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

% u_vecall = [0; 0.1; 0.2; 0.3; 0.5; 0.8];
% length_variable = length(u_vecall);

% statistical variable initialize
Falsemat_tau = zeros(MC, length_variable);
Falsemat_CA = zeros(MC, length_variable);

Detectmat_tau = zeros(MC, length_variable);
Detectmat_CA = zeros(MC, length_variable);

% Mento Carlo method
tic;

for sp_idx = 1 : length_variable
    SNR = SNRvec_all(sp_idx);
    % u_set = 0;

    for mc_idx = 1 : MC
        handle_waitbar = waitbar(((sp_idx - 1) * MC + mc_idx) / ...
            (MC * length_variable));

        omega_true = zeros(K, 1);
        omega_min = 2 * pi / Nx;
        sigma_vec = sigma_n * exprnd(mu_sigma, [Nx, S_snap]);
        y_noise = sqrt(sigma_vec / 2) .* (randn(Nx, S_snap) + ...
            1j * randn(Nx, S_snap));
        omega_true(1) = pi * (2 * rand - 1);
        for k = 2 : K
            th = pi * (2 * rand - 1);
            while min(abs((wrapToPi(th - omega_true(1 : k - 1))))) < ...
                omega_min
                th = pi * (2 * rand - 1);
            end
            omega_true(k) = th;
        end
        omega_true = wrapTo2Pi(omega_true);
        gain_true = bsxfun(@times, sqrt(sigma_n) * (10 .^ (SNR / 20)), ...
            exp(1j*2*pi*rand(K, S_snap)));  % K * S_snap
        y_full = exp(1j * (0 : (Nx - 1)).' * omega_true.') / ...
            sqrt(Nx) * gain_true;
        y = Smat_com * y_full + y_noise;

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
    end
end


toc;
delete(handle_waitbar);

% after care
Falserate_tau = mean(Falsemat_tau);
Detectrate_tau = mean(Detectmat_tau);

Falserate_CA = mean(Falsemat_CA);
Detectrate_CA = mean(Detectmat_CA);

if MC > 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PFAvsSigmavar.mat'];
    save(filename_now, 'Nx', 'P_oe', 'K', 'u_vecall', 'length_variable',...
    'Falsemat_tau', 'Falsemat_CA', 'Detectmat_tau', 'Detectmat_CA');
end

% plot the result
lw = 2;
fsz = 12;
msz = 8;


figure(1)
plot(SNRvec_all, P_oe * ones(1, length_variable), '--k', 'Linewidth', lw)
hold on;
plot(SNRvec_all, Falserate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
plot(SNRvec_all, Falserate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
legend('$\bar{\rm P}_{\rm FA} = 0.01$', 'NOMP', ...
    'NOMP-CFAR', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('Strength of Noise fluctuation u (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured $\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(2)
plot(SNRvec_all, Detectrate_tau, '-ro', 'Linewidth', lw, 'Markersize', msz)
hold on;
plot(SNRvec_all, Detectrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
xlabel('Strength of Noise Fluctuation u (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('Measured ${\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)
