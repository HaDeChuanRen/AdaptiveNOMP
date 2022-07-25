% last update: 2022/7/24
% test the performance of NMOP-CFAR-2D

clc; clear; close all;

rng(5)
MC = 100;

% Define Scenario
Nx = 64; % Length of Sinusoid
My = 32; % width of Sinusoid
NM_num = Nx * My;
K = 8;
sigma_n = 1;              % noise variance sigma^2, instead of sigma

S_snap = 1;
SNRvec_all = 6 : 2 : 26;
length_SNR = length(SNRvec_all);
array_Fun = @(omega, N) exp(1j * (0 : (N - 1))' * omega) / sqrt(N);


% algorithm parameters
K_max = K * 2;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;

guard_n = 2;
guard_m = 3;
training_n = 5;
training_m = 5;
guard_training_size = [guard_n, guard_m, training_n, training_m];
N_r = (2 * training_n + 1) * (2 * training_m + 1) - (2 * guard_n + 1) * (2 * guard_m + 1);

alpha_set = alpha_PoebyS(P_oe, NM_num, N_r, S_snap);

% statistical variable initialize
Overestmat_CA = zeros(MC, length_SNR);
Equalmat_CA = zeros(MC, length_SNR);

% Mento-Carlo method
tic;

for sp_idx = 1 : length_SNR
    SNR = SNRvec_all(sp_idx);

    for mc_idx = 1 : MC
        handle_waitbar = waitbar(((sp_idx - 1) * MC + mc_idx) / (MC * length_SNR));
        omega_true = (2 * rand(K, 2) - 1) * pi;
        gain_phi = rand(K, 1) * 2 * pi;
        gain_true = gain_phi * sqrt(sigma_n) * (10 .^ (SNR / 20));

        y_vec = zeros(NM_num, 1);
        for k_idx = 1 : K
            y_veck = kron(array_Fun(omega_true(k_idx, 2), My),...
            array_Fun(omega_true(k_idx, 1), Nx));
            y_vec = gain_true(k_idx) * y_veck + y_vec;
        end
        y_matrix_pure = reshape(y_vec, Nx, My);
        ymat= y_matrix_pure + (sigma_n / 2) * (1j * randn(Nx, My) + randn(Nx, My));

        [omegavec_CA, gainvec_CA, ~] = ...
        NOMP2D_CFAR(ymat, alpha_set, guard_training_size, K_max);
        if length(gainvec_CA) > K
            Overestmat_CA(mc_idx, sp_idx) = 1;
        elseif length(gainvec_CA) == K
            Equalmat_CA(mc_idx, sp_idx) = 1;
        end

    end
end

toc;
delete(handle_waitbar);

% after care
Overestrate_CA = mean(Overestmat_CA);
Equalrate_CA = mean(Equalmat_CA);

% plot the result
lw = 1.6;
fsz = 12;
msz = 10;


figure(1)
plot(SNRvec_all, P_oe * ones(1, length_SNR), '--k', 'Linewidth', lw)
hold on;
plot(SNRvec_all, Overestrate_CA, '-b+', 'Linewidth', lw, 'Markersize', msz)
legend('${\rm P}_{\rm OE} = 0.01$', 'NOMP-CA', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('${\rm SNR}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(2)
plot(SNRvec_all, Equalrate_CA, '-b+', 'Linewidth', lw)
xlabel('${\rm SNR}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured ${\rm P}_{\rm equal}$', 'Interpreter', 'latex', 'Fontsize', fsz)
