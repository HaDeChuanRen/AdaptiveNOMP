% last updated date: 2022.7.13
% 2022.6.25: 1. compare the Poe performance of NOMP and CA-NOMP

clc; clear; close all;
set(0,'DefaultLineMarkerSize',12);
set(0,'DefaultTextFontSize',14);
set(0,'DefaultAxesFontSize',14);


rng(5);
MC = 5000;

% Define Scenario
N = 256; % Length of Sinusoid
% number of sinusoids in the mixture of sinusoids
% K = 8;
K_max_all = [1, 4, 8];


% SNR = 22;
% SNR_all = [20, 22, 24];
CFAR_method = 'CA';
gamma_oversamping = 4;

% guard_size = 5;
% training_size = 30;
% guard_training_size = [guard_size, training_size];
% training_size * 2
N_r = 60;

% training_size10 = 10;
% guard_training_size10 = [guard_size, training_size10];

sigma_n = 1;              % noise variance sigma^2, instead of sigma
T = 1;


Poe_all = [1e-3, 2e-3, 3e-3, 5e-3, 7e-3, 1e-2, 2e-2, 3e-2, 5e-2, 7e-2, 0.1];

% Poe_all = [1e-2, 2e-2, 3e-2, 5e-2, 0.1];
% cauculate the alpha of the NOMP-CA
N_alpha_g = 10000;
P_OE_bar = 0.01;
alpha_true = -log(1-(1-P_OE_bar)^(1/N));
alpha_grid = linspace(alpha_true/30,5*alpha_true,N_alpha_g);

res = zeros(N_alpha_g,1);
for grid_alpha = 1:N_alpha_g
    fun1 = @(xvar) exp(N*log(1-exp(-alpha_grid(grid_alpha)/(2*N_r)*xvar))+...
    (N_r-1)*log(xvar)-xvar/2-N_r*log(2)-sum(log(1:N_r-1)));
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


% omega_true = zeros(K, 1);
% omega_min = 2 * pi / N;
% normal measurements
% M = N; % number of measurements = N
% S = eye(N);

M = round(N / 2);
% type of measurement matrix
measure_type = 'cmplx_bernoulli';
% windowing weights
% options 'all_ones' (default), 'hamming' and 'hann'
window_type = 'all_ones';
S = generateMeasMat(N, M, measure_type, window_type);

Overest_tau = zeros(MC, Num_Poe);
Overest_CA_small = zeros(MC, Num_Poe);
Overest_CA_medium = zeros(MC, Num_Poe);
Overest_CA_big = zeros(MC, Num_Poe);



tic;
for sp_idx = 1 : Num_Poe
    P_oe = Poe_all(sp_idx);

    alpha_set = alpha_CA(sp_idx);
    tau_NOMP = sigma_n * chi2inv((1 - P_oe) ^ (1 / N), 2 * T) / 2;

    for mc = 1 : MC
        handle_waitbar = waitbar(((sp_idx - 1) * MC + mc) / (MC * Num_Poe));

        y_noise = sqrt(sigma_n / 2) * (randn(M, T) + 1j*randn(M, T));

        [omegaList_tau, gainList_tau, ~] = MNOMP(y_noise, S, tau_NOMP);
        if ~isempty(omegaList_tau)
            Overest_tau(mc, sp_idx) = 1;
        end

        [omegaList_CA_small, gainList_CA_small, ~] =...
        MNOMP_CFAR_alpha(y_noise, S, alpha_set, N_r, K_max_all(1), CFAR_method);
        if ~isempty(omegaList_CA_small)
            Overest_CA_small(mc, sp_idx) = 1;
        end

        % [omegaList_CA_medium, gainList_CA_medium, ~] =...
        % MNOMP_CFAR_alpha(y_noise, S, alpha_set, N_r, K_max_all(2), CFAR_method);
        % if ~isempty(omegaList_CA_medium)
        %     Overest_CA_medium(mc, sp_idx) = 1;
        % end

        % [omegaList_CA_big, gainList_CA_big, ~] =...
        % MNOMP_CFAR_alpha(y_noise, S, alpha_set, N_r, K_max_all(3), CFAR_method);
        % if ~isempty(omegaList_CA_big)
        %     Overest_CA_big(mc, sp_idx) = 1;
        % end

    end
end

toc;
delete(handle_waitbar);
Overest_rate_tau = mean(Overest_tau);
Overest_rate_CA_small = mean(Overest_CA_small);
Overest_rate_CA_medium = mean(Overest_CA_medium);
Overest_rate_CA_big = mean(Overest_CA_big);




lw = 1.6;
fsz = 12;
msz = 10;

figure(1);
loglog(Poe_all, Poe_all, '--k','Linewidth',lw)
hold on;
loglog(Poe_all, Overest_rate_tau,'ro','Linewidth',lw,'Markersize',msz)
loglog(Poe_all, Overest_rate_CA_small,'b+','Linewidth',lw,'Markersize',msz)
% loglog(Poe_all, Overest_rate_CA_medium,'b^','Linewidth',lw,'Markersize',msz)
% loglog(Poe_all, Overest_rate_CA_big, 'bv', 'Linewidth',lw,'Markersize',msz)
legend('${\rm P}_{\rm OE}(nominal)$', 'NOMP', 'CA-NOMP($K_{\rm max}$ = 1)',...
'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('nominal ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)

% 'CA-NOMP($K_{\rm max}$ = 4)', 'CA-NOMP($K_{\rm max}$ = 8)',...

if MC >= 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PoevsSNR.mat'];
    save(filename_now, 'N', 'Poe_all', 'Overest_tau', 'Overest_CA_small',...
    'Overest_CA_medium', 'Overest_CA_big');
end

