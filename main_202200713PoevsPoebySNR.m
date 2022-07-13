% last updated date: 2022.7.13
% 2022.6.25: 1. compare the Poe performance of NOMP and CA-NOMP

clc; clear; close all;
set(0,'DefaultLineMarkerSize',12);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',14);


rng(5);
MC = 5;

% Define Scenario
N = 256; % Length of Sinusoid
% number of sinusoids in the mixture of sinusoids
K = 16;
K_max = K * 2;


% SNR = 22;
SNR_all = [20, 22, 24];
CFAR_method = 'CA';
gamma_oversamping = 4;

guard_size = 5;
training_size = 30;
guard_training_size = [guard_size, training_size];
N_r = training_size * 2;

% training_size10 = 10;
% guard_training_size10 = [guard_size, training_size10];

sigma_n = 1;              % noise variance sigma^2, instead of sigma
T = 1;


% Poe_all = [1e-3, 2e-3, 3e-3, 5e-3, 1e-2, 2e-2, 3e-2, 5e-2, 0.1];

Poe_all = [1e-2, 2e-2, 3e-2, 5e-2, 0.1];
% cauculate the alpha of the NOMP-CA
N_alpha_g = 10000;
P_OE_bar = 0.01;
alpha_true = -log(1-(1-P_OE_bar)^(1/N));
alpha_grid = linspace(alpha_true/30,5*alpha_true,N_alpha_g);
Nr = 50;
res = zeros(N_alpha_g,1);
for grid_alpha = 1:N_alpha_g
    fun1 = @(xvar) exp(N*log(1-exp(-alpha_grid(grid_alpha)/(2*Nr)*xvar))+(Nr-1)*log(xvar)-xvar/2-Nr*log(2)-sum(log(1:Nr-1)));
    res(grid_alpha) = integral(fun1,0,Inf);
end
Num_Poe = length(Poe_all);
alpha_CA = zeros(length(Poe_all),1);
% nominal_POE = 0.01:0.01:0.12;
oneminusPOE = 1-Poe_all;
for idx = 1:length(Poe_all)
    oneminusPOE_idx = oneminusPOE(idx);
    [~,idxmin] = min(abs(oneminusPOE_idx-res));
    alpha_CA(idx) = alpha_grid(idxmin);
end


omega_true = zeros(K, 1);
omega_min = 2 * pi / N;
% normal measurements
M = N; % number of measurements = N
S = eye(N);

% Threshold_collect = zeros(MC, Num_Poe);
% tau_NOMP_collect = zeros(1, Num_Poe);

Overest_tau_SNRlow = zeros(MC, Num_Poe);
% Overest_tau_SNRmedium = zeros(MC, Num_Poe);
% Overest_tau_SNRhigh = zeros(MC, Num_Poe);

Overest_CA_SNRlow = zeros(MC, Num_Poe);
% Overest_CA_SNRmedium = zeros(MC, Num_Poe);
% Overest_CA_SNRhigh = zeros(MC, Num_Poe);


tic;
for sp_idx = 1 : Num_Poe
    P_oe = Poe_all(sp_idx);

    alpha_set = alpha_CA(sp_idx);
    tau_NOMP = sigma_n * chi2inv((1 - P_oe) ^ (1 / N), 2 * T) / 2;

    for mc = 1 : MC
        handle_waitbar = waitbar(((sp_idx - 1) * MC + mc) / (MC * Num_Poe));
        % SNR = SNR_min + SNR_delta * rand(K, 1);

        % [y, omega_true, gain_true] = create_yvector(K, T, N, SNR, sigma_n, S);
        omega_true = zeros(K, 1);
        omega_true(1) = pi * (2 * rand - 1);
        for k = 2 : K
            th = pi * (2 * rand - 1);
            while min(abs((wrapToPi(th - omega_true(1 : k - 1))))) < omega_min
                th = pi * (2*rand - 1);
            end
            omega_true(k) = th;
        end
        omega_true = wrapTo2Pi(omega_true);
        gain_SNRlow = bsxfun(@times, sqrt(sigma_n) * (10.^(SNR_all(1) / 20)),...
        exp(1j*2*pi*rand(K,T)));  % K*T
        % gain_SNRmedium = bsxfun(@times, sqrt(sigma_n) * (10.^(SNR_all(2) / 20)),...
        % exp(1j*2*pi*rand(K,T)));  % K*T
        % gain_SNRhigh = bsxfun(@times, sqrt(sigma_n) * (10.^(SNR_all(3) / 20)),...
        % exp(1j*2*pi*rand(K,T)));  % K*T
        % original signal
        % y_full = exp(1j * (0:(N-1)).' * omega_true.')/sqrt(N) * gain_true;
        % y_full = zeros(N, T);
        noise = sqrt(sigma_n / 2) * (randn(N, T) + 1j*randn(N, T));
        % y_noisy = S * y_full + noise;
        y_SNRlow = exp(1j * (0:(N-1)).' * omega_true.')/sqrt(N) * gain_SNRlow + noise;
        % y_SNRmedium = exp(1j * (0:(N-1)).' * omega_true.')/sqrt(N) * gain_SNRmedium + noise;
        % y_SNRhigh = exp(1j * (0:(N-1)).' * omega_true.')/sqrt(N) * gain_SNRhigh + noise;

        [omegaList_tau_SNRlow, ~, ~] = MNOMP(y_SNRlow, S, tau_NOMP);
        if length(omegaList_tau_SNRlow) > K
            Overest_tau_SNRlow(mc, sp_idx) = 1;
        end

        % [omegaList_tau_SNRmedium, ~, ~] =...
        % MNOMP_CFAR_alpha(y_SNRmedium, S, alpha_set, guard_training_size, K_max, CFAR_method);
        % if length(omegaList_tau_SNRmedium) > K
        %     Overest_tau_SNRmedium(mc, sp_idx) = 1;
        % end

        % [omegaList_tau_SNRhigh, ~, ~] =...
        % MNOMP_CFAR_alpha(y_SNRhigh, S, alpha_set, guard_training_size, K_max, CFAR_method);
        % if length(omegaList_tau_SNRhigh) > K
        %     Overest_tau_SNRhigh(mc, sp_idx) = 1;
        % end

        [omegaList_CA_SNRlow, ~, ~] =...
        MNOMP_CFAR_alpha(y_SNRlow, S, alpha_set, N_r, K_max, CFAR_method) ;
        if length(omegaList_CA_SNRlow) > K
            Overest_CA_SNRlow(mc, sp_idx) = 1;
        end

        % [omegaList_CA_SNRmedium, ~, ~] =...
        % MNOMP_CFAR_alpha(y_SNRmedium, S, alpha_set, N_r, K_max, CFAR_method);
        % if length(omegaList_CA_SNRmedium) > K
        %     Overest_CA_SNRmedium(mc, sp_idx) = 1;
        % end

        % [omegaList_CA_SNRhigh, ~, ~] =...
        % MNOMP_CFAR_alpha(y_SNRhigh, S, alpha_set, guard_training_size, K_max, CFAR_method);
        % if length(omegaList_CA_SNRhigh) > K
        %     Overest_CA_SNRhigh(mc, sp_idx) = 1;
        % end

    end
end

toc;
delete(handle_waitbar);
Overest_rate_tau_SNRlow = mean(Overest_tau_SNRlow);
% Overest_rate_tau_SNRmedium = mean(Overest_tau_SNRmedium);
% Overest_rate_tau_SNRhigh = mean(Overest_tau_SNRhigh);

Overest_rate_CA_SNRlow = mean(Overest_CA_SNRlow);
% Overest_rate_CA_SNRmedium = mean(Overest_CA_SNRmedium);
% Overest_rate_CA_SNRhigh = mean(Overest_CA_SNRhigh);



lw = 1.6;
fsz = 12;
msz = 10;
figure(1);
loglog(Poe_all, Poe_all, '--k','Linewidth',lw)
hold on;
loglog(Poe_all, Overest_rate_tau_SNRlow,'ro','Linewidth',lw,'Markersize',msz)
loglog(Poe_all, Overest_rate_CA_SNRlow,'bo','Linewidth',lw,'Markersize',msz)
% loglog(Poe_all, Overest_rate_tau_SNRmedium,'r+','Linewidth',lw,'Markersize',msz)
% loglog(Poe_all, Overest_rate_CA_SNRmedium,'b+','Linewidth',lw,'Markersize',msz)
% loglog(Poe_all, Overest_rate_tau_SNRhigh,'r^','Linewidth',lw,'Markersize',msz)
% loglog(Poe_all, Overest_rate_CA_SNRhigh,'b^','Linewidth',lw,'Markersize',msz)
legend('${\rm P}_{\rm OE} = 0.01$', 'SNR = 20dB','Interpreter', 'latex', 'Fontsize', fsz)
xlabel('nominal ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)
% title_SNR = ['SNR=', num2str(SNR), ', N_r=', num2str(N_r)];
% title(title_SNR)

% , 'SNR = 16dB',...
% 'SNR = 18dB', 

filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDvsSNR.mat'];
save(filename_now, 'N', 'Poe_all', 'Overest_rate_CA_SNRlow', 'Overest_rate_tau_SNRlow');

% 'False_rate_infty', 'Overest_rate_infty', 'Detect_rate_infty',...
% 'Equal_rate_infty', 'Miss_rate_infty', 'MSEdB_infty'
% 'Overest_rate_CA_SNRmedium', 'Overest_rate_CA_SNRhigh'
