% last update: 2022.7.7
% 2022.7.6 : set the omega_min as 2.5 * 2 * pi / N;
% 2022.7.3 : change the save data as MC * Num_SNR

clc; clear; close all;
% set(0,'DefaultLineMarkerSize',12);
% set(0,'DefaultTextFontSize',14);
% set(0,'DefaultAxesFontSize',14);


rng(5);
MC = 3000;

% Define Scenario
N = 256; % Length of Sinusoid

K = 8;
K_max = 2 * K;
SNR_min_all = 12 : 1 : 22;
Num_SNR = length(SNR_min_all);
SNR_delta = 0;
guard_size = 8;
training_size = 50;

sigma_n = 1;              % noise variance sigma^2, instead of sigma
S_snap = 1;
P_oe001 = 0.01;
P_oe005 = 0.05;

N_r50 = 50;
N_r20 = 20;
CFAR_method = 'CA';

gamma_oversamping = 4;

tau_NOMP_001Poe = sigma_n * chi2inv((1 - P_oe001) ^ (1 / N), 2 * S_snap) / 2;
alpha_CA25_001Poe = alpha_PoebyS(P_oe001, N, N_r50, S_snap);


omega_true = zeros(K, 1);
omega_min = 2.5 * 2 * pi / N;
% normal measurements
M = N; % number of measurements = N
S_matrix = eye(N);

% M = round(N / 2);
% % type of measurement matrix
% measure_type = 'cmplx_bernoulli';
% % windowing weights
% % options 'all_ones' (default), 'hamming' and 'hann'
% window_type = 'all_ones';
% S = generateMeasMat(N, M, measure_type, window_type);

Miss_tau_001Poe = zeros(MC, Num_SNR);
False_tau_001Poe = zeros(MC, Num_SNR);
Overest_tau_001Poe = zeros(MC, Num_SNR);
Detect_tau_001Poe = zeros(MC, Num_SNR);
Equal_tau_001Poe = zeros(MC, Num_SNR);
Error_ave_tau_001Poe = zeros(MC, Num_SNR);
reconErr_tau_001Poe = zeros(MC, Num_SNR);


Miss_CA25_001Poe = zeros(MC, Num_SNR);
False_CA25_001Poe = zeros(MC, Num_SNR);
Overest_CA25_001Poe = zeros(MC, Num_SNR);
Detect_CA25_001Poe = zeros(MC, Num_SNR);
Equal_CA25_001Poe = zeros(MC, Num_SNR);
Error_ave_CA25_001Poe = zeros(MC, Num_SNR);
reconErr_CA25_001Poe= zeros(MC, Num_SNR);

Miss_VALSE = zeros(MC, Num_SNR);
False_VALSE = zeros(MC, Num_SNR);
Overest_VALSE = zeros(MC, Num_SNR);
Detect_VALSE = zeros(MC, Num_SNR);
Equal_VALSE = zeros(MC, Num_SNR);
Error_ave_VALSE = zeros(MC, Num_SNR);
reconErr_VALSE= zeros(MC, Num_SNR);



tic;
for mc = 1 : MC
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
    noise = sqrt(sigma_n / 2) * (randn(M, S_snap) + 1j*randn(M, S_snap));
    gain_phi = exp(1j * 2 * pi * rand(K, S_snap));

    for sp_idx = 1 : Num_SNR
        SNR_min = SNR_min_all(sp_idx);
        handle_waitbar = waitbar(((mc - 1) * Num_SNR + sp_idx) / (MC * Num_SNR));
        SNR = SNR_min + SNR_delta * rand(K, 1);

        % [y, omega_true, gain_true] = create_yvector(K, T, N, SNR, sigma_n, S);
        gain_true = bsxfun(@times, sqrt(sigma_n) * (10.^(SNR / 20)), gain_phi);  % K*T
        % original signal
        y_full = exp(1j * (0:(N-1)).' * omega_true.')/sqrt(N) * gain_true;
        y_noisy = S_matrix * y_full + noise;
        y = y_noisy;

        % if (mc ~= 434) || (sp_idx ~= 9)
        %     continue;
        % end

        % if (mc ~= 37) || (sp_idx ~= 1)
        %     continue;
        % end

        % NOMP P_oe = 0.01
        [omegaList_tau_001Poe, gainList_tau_001Poe, ~] = MNOMP(y, S_matrix, tau_NOMP_001Poe);
        results_struct_tau_001Poe = False_Detection(omega_true, gain_true,...
        omegaList_tau_001Poe, gainList_tau_001Poe, N);
        False_tau_001Poe(mc, sp_idx) = results_struct_tau_001Poe.False_Eve;
        Overest_tau_001Poe(mc, sp_idx) = results_struct_tau_001Poe.Overest_Eve;
        Detect_tau_001Poe(mc, sp_idx) = results_struct_tau_001Poe.Detect_Eve;
        Equal_tau_001Poe(mc, sp_idx) = results_struct_tau_001Poe.Equal_Eve;
        Error_ave_tau_001Poe(mc, sp_idx) = results_struct_tau_001Poe.error_ave;
        Miss_tau_001Poe(mc, sp_idx) = results_struct_tau_001Poe.Miss_Eve;
        reconErr_tau_001Poe(mc, sp_idx) = results_struct_tau_001Poe.reconErr;

        % NOMP-CA P_oe = 0.01 N_r = 50
        [omegaList_CA25_001Poe, gainList_CA25_001Poe, ~] =...
        MNOMP_CFAR_alpha(y, S_matrix, alpha_CA25_001Poe, training_size, K_max, CFAR_method);
        results_struct_CA25_001Poe = False_Detection(omega_true, gain_true,...
        omegaList_CA25_001Poe, gainList_CA25_001Poe, N);
        False_CA25_001Poe(mc, sp_idx) = results_struct_CA25_001Poe.False_Eve;
        Overest_CA25_001Poe(mc, sp_idx) = results_struct_CA25_001Poe.Overest_Eve;
        Detect_CA25_001Poe(mc, sp_idx) = results_struct_CA25_001Poe.Detect_Eve;
        Equal_CA25_001Poe(mc, sp_idx) = results_struct_CA25_001Poe.Equal_Eve;
        Error_ave_CA25_001Poe(mc, sp_idx) = results_struct_CA25_001Poe.error_ave;
        Miss_CA25_001Poe(mc, sp_idx) = results_struct_CA25_001Poe.Miss_Eve;
        reconErr_CA25_001Poe(mc, sp_idx) = results_struct_CA25_001Poe.reconErr;

        out_VALSE = MVALSE_best(y, (0 : (N-1))', 2, y_full);
        omegaList_VALSE = wrapTo2Pi(out_VALSE.freqs);
        gainList_VALSE = out_VALSE.amps;
        results_struct_VALSE = False_Detection(omega_true, gain_true,...
        omegaList_VALSE, gainList_VALSE, N);
        False_VALSE(mc, sp_idx) = results_struct_VALSE.False_Eve;
        Overest_VALSE(mc, sp_idx) = results_struct_VALSE.Overest_Eve;
        Detect_VALSE(mc, sp_idx) = results_struct_VALSE.Detect_Eve;
        Equal_VALSE(mc, sp_idx) = results_struct_VALSE.Equal_Eve;
        Error_ave_VALSE(mc, sp_idx) = results_struct_VALSE.error_ave;
        Miss_VALSE(mc, sp_idx) = results_struct_VALSE.Miss_Eve;
        reconErr_VALSE(mc, sp_idx) = out_VALSE.mse(end);
    end
end

toc;
delete(handle_waitbar);
SNR_delta = 0;
CRB_omega = 10 * log10((6 / (N^2 - 1)) * 10 .^ (- (SNR_min_all + SNR_delta / 2) / 10));
Num_SNR = length(SNR_min_all);

False_rate_tau_001Poe = mean(False_tau_001Poe);
Overest_rate_tau_001Poe = mean(Overest_tau_001Poe);
Detect_rate_tau_001Poe = mean(Detect_tau_001Poe);
Equal_rate_tau_001Poe = mean(Equal_tau_001Poe);
Miss_rate_tau_001Poe = mean(Miss_tau_001Poe);
MSEdB_tau_001Poe = 10 * log10(sum(Error_ave_tau_001Poe) ./ sum(Equal_tau_001Poe));
MSEdB_tau_001Poe(Equal_rate_tau_001Poe < 0.3) = nan;
reconErr_tau_001Poe_ave = 10 * log10(mean(reconErr_tau_001Poe));

False_rate_CA25_001Poe = mean(False_CA25_001Poe);
Overest_rate_CA25_001Poe = mean(Overest_CA25_001Poe);
Detect_rate_CA25_001Poe = mean(Detect_CA25_001Poe);
Equal_rate_CA25_001Poe = mean(Equal_CA25_001Poe);
Miss_rate_CA25_001Poe = mean(Miss_CA25_001Poe);
MSEdB_CA25_001Poe = 10 * log10(sum(Error_ave_CA25_001Poe) ./ sum(Equal_CA25_001Poe));
MSEdB_CA25_001Poe(Equal_rate_CA25_001Poe < 0.3) = nan;
reconErr_CA25_001Poe_ave = 10 * log10(mean(reconErr_CA25_001Poe));

False_rate_VALSE = mean(False_VALSE);
Overest_rate_VALSE = mean(Overest_VALSE);
Detect_rate_VALSE = mean(Detect_VALSE);
Equal_rate_VALSE = mean(Equal_VALSE);
Miss_rate_VALSE = mean(Miss_VALSE);
MSEdB_VALSE = 10 * log10(sum(Error_ave_VALSE) ./ sum(Equal_VALSE));
MSEdB_VALSE(Equal_rate_VALSE < 0.3) = nan;
reconErr_VALSE_ave = 10 * log10(mean(reconErr_VALSE));


lw = 2;
fsz = 12;
msz = 8;

% figure(1);
% plot(SNR_min_all, P_oe001 * ones(1, Num_SNR), '--k','Linewidth',lw)
% hold on;
% plot(SNR_min_all, Overest_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
% plot(SNR_min_all, Overest_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
% legend('${\rm P}_{\rm OE} = 0.01$', 'NOMP', ...
%     'NOMP-CFAR $(N_r = 50)$', 'Interpreter', 'latex', 'Fontsize', fsz)
% xlabel('${\rm SNR}$ (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
% ylabel('measured ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)

figure(1);
plot(SNR_min_all, P_oe001 * ones(1, Num_SNR), '--k','Linewidth',lw)
hold on;
plot(SNR_min_all, False_rate_VALSE,'-md','Linewidth', lw,'Markersize', msz)
plot(SNR_min_all, False_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, False_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
legend('$\bar{\rm P}_{\rm FA} = 0.01$', 'VALSE', 'NOMP', ...
    'NOMP-CFAR', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('Measured $\bar{\rm P}_{\rm FA}$','Interpreter','latex','Fontsize',fsz)

figure(4);
semilogy(SNR_min_all, P_oe001 * ones(1, Num_SNR), '--k','Linewidth',lw)
hold on;
semilogy(SNR_min_all, False_rate_VALSE,'-md','Linewidth', lw,'Markersize', msz)
semilogy(SNR_min_all, False_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
semilogy(SNR_min_all, False_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
legend('$\bar{\rm P}_{\rm FA} = 0.01$', 'VALSE', 'NOMP', ...
    'NOMP-CFAR', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('Measured $\bar{\rm P}_{\rm FA}$','Interpreter','latex','Fontsize',fsz)


figure(2);
hold on;
plot(SNR_min_all, Equal_rate_VALSE,'-md','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, Equal_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, Equal_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('Measured ${\rm P}(\hat{K} = K)$','Interpreter','latex','Fontsize',fsz)

figure(3);
hold on;
plot(SNR_min_all, Detect_rate_VALSE,'-md','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, Detect_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, Detect_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('Measured $\bar{\rm P}_{\rm D}$','Interpreter','latex','Fontsize',fsz)



figure(5);
plot(SNR_min_all, CRB_omega,'--k','Linewidth',lw)
hold on;
plot(SNR_min_all, MSEdB_VALSE,'-md','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, MSEdB_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, MSEdB_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('MSE($\hat{\mathbf \omega}$) (dB)','Interpreter','latex','Fontsize',fsz)
legend('CRB','Fontsize',fsz)

% figure(6);
% plot(SNR_min_all, Miss_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
% hold on;
% plot(SNR_min_all, Miss_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
% xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
% ylabel('measured ${\rm P}_{\rm M}$','Interpreter','latex','Fontsize',fsz)

figure(7);
hold on;
plot(SNR_min_all, reconErr_VALSE_ave,'-md','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, reconErr_tau_001Poe_ave,'-ro','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, reconErr_CA25_001Poe_ave,'-b+','Linewidth',lw,'Markersize',msz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('NMSE($\hat{\mathbf{z}}$) (dB)','Interpreter','latex','Fontsize',fsz)
% \|\hat{\mathbf{z}} - \mathbf{z}\|

if MC > 100
    filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PDPFAMSEvsSNR.mat'];
    save(filename_now, 'N', 'P_oe001', 'K', 'SNR_min_all',  'CRB_omega',...
    'False_tau_001Poe', 'Overest_tau_001Poe', 'Detect_tau_001Poe',...
    'Equal_tau_001Poe', 'Miss_tau_001Poe', 'Error_ave_tau_001Poe', 'reconErr_tau_001Poe',...
    'False_CA25_001Poe', 'Overest_CA25_001Poe', 'Detect_CA25_001Poe',...
    'Equal_CA25_001Poe', 'Miss_CA25_001Poe', 'Error_ave_CA25_001Poe', 'reconErr_CA25_001Poe', ...
    'False_VALSE', 'Overest_VALSE', 'Detect_VALSE',...
    'Equal_VALSE', 'Miss_VALSE', 'Error_ave_VALSE', 'reconErr_VALSE');
end





