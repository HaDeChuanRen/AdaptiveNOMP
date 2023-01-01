clc; clear; close all;

N = 256; % Length of Sinusoid
% number of sinusoids in the mixture of sinusoids
K = 16;
S_snap = 1;
SNR_all = [14, 15, 18];
sigma_n = 1;
gain_low = sqrt(sigma_n) * (10 ^ (SNR_all(1) / 20));
gain_medium = sqrt(sigma_n) * (10 ^ (SNR_all(2) / 20));
gain_high = sqrt(sigma_n) * (10 ^ (SNR_all(3) / 20));
N_r = 60;

beta_ave = 0.88;

P_FA_all = 0.01 : 0.01 : 0.12;
Num_Poe = length(P_FA_all);

PDveclow_tau = zeros(Num_Poe, 1);
PDvecmedium_tau = zeros(Num_Poe, 1);
PDvechigh_tau = zeros(Num_Poe, 1);
PDveclow_CFAR = zeros(Num_Poe, 1);
PDvecmedium_CFAR = zeros(Num_Poe, 1);
PDvechigh_CFAR = zeros(Num_Poe, 1);

for sp_idx = 1 : Num_Poe
    P_oe = P_FA_all(sp_idx);
    alpha_set = alpha_PoebyS(P_oe, N, N_r);
    tau_NOMP = sigma_n * chi2inv((1 - P_oe) ^ (1 / N), 2 * S_snap) / 2;

    PD_low = marcumq(sqrt(2 / sigma_n) * gain_low * beta_ave, ...
    sqrt(2 * tau_NOMP / sigma_n));
    PDveclow_tau(sp_idx) = PD_low ^ K;
    fun_PDlow = @(Tvar) marcumq(sqrt(2 / sigma_n) * gain_low * beta_ave, ...
        sqrt(2 * Tvar)) .* exp((N_r - 1) * log(Tvar) - (N_r / alpha_set) * Tvar - ...
    sum(log(1 : N_r - 1)) + N_r * log(N_r / alpha_set));
    PDveclow_CFAR(sp_idx) = integral(fun_PDlow, 0, Inf) ^ K;

    PD_medium = marcumq(sqrt(2 / sigma_n) * gain_medium * beta_ave, ...
    sqrt(2 * tau_NOMP / sigma_n));
    PDvecmedium_tau(sp_idx) = PD_medium ^ K;
    fun_PDmedium = @(Tvar) marcumq(sqrt(2) * gain_medium * beta_ave / sigma_n, ...
        sqrt(2 * Tvar)) .* exp((N_r - 1) * log(Tvar) - (N_r / alpha_set) * Tvar - ...
    sum(log(1 : N_r - 1)) + N_r * log(N_r / alpha_set));
    PDvecmedium_CFAR(sp_idx) = integral(fun_PDmedium, 0, Inf) ^ K;

    PD_high = marcumq(sqrt(2 / sigma_n) * gain_high * beta_ave, ...
    sqrt(2 * tau_NOMP / sigma_n));
    PDvechigh_tau(sp_idx) = PD_high ^ K;
    fun_PDhigh = @(Tvar) marcumq(sqrt(2 / sigma_n) * gain_high * beta_ave, ...
        sqrt(2 * Tvar)) .* exp((N_r - 1) * log(Tvar) - (N_r / alpha_set) * Tvar - ...
    sum(log(1 : N_r - 1)) + N_r * log(N_r / alpha_set));
    PDvechigh_CFAR(sp_idx) = integral(fun_PDhigh, 0, Inf) ^ K;

end

% T_var = 0 : 1000;
% figure;
% plot(T_var, fun_PDlow(T_var));

load('20220909T024133_mc3000_PfavsPfabySNR.mat');
False_rate_tau_SNRlow = mean(False_tau_SNRlow);
Detect_rate_tau_SNRlow = mean(Detect_tau_SNRlow);
False_rate_tau_SNRmedium = mean(False_tau_SNRmedium);
Detect_rate_tau_SNRmedium = mean(Detect_tau_SNRmedium);
False_rate_tau_SNRhigh = mean(False_tau_SNRhigh);
Detect_rate_tau_SNRhigh = mean(Detect_tau_SNRhigh);

False_rate_CA_SNRlow = mean(False_CA_SNRlow);
Detect_rate_CA_SNRlow = mean(Detect_CA_SNRlow);
False_rate_CA_SNRmedium = mean(False_CA_SNRmedium);
Detect_rate_CA_SNRmedium = mean(Detect_CA_SNRmedium);
False_rate_CA_SNRhigh = mean(False_CA_SNRhigh);
Detect_rate_CA_SNRhigh = mean(Detect_CA_SNRhigh);

lw = 2;
fsz = 12;
msz = 8;

figure;
plot(Poe_all * 100, Poe_all * 100, '--k', 'Linewidth', lw)
hold on;
plot(Poe_all * 100, False_rate_tau_SNRlow * 100, 'ro', 'Linewidth', lw,'Markersize',msz)
plot(Poe_all * 100, False_rate_CA_SNRlow * 100, 'bo','Linewidth',lw,'Markersize',msz)
plot(Poe_all * 100, False_rate_tau_SNRmedium * 100,'r+','Linewidth',lw,'Markersize',msz)
plot(Poe_all * 100, False_rate_CA_SNRmedium * 100,'b+','Linewidth',lw,'Markersize',msz)
plot(Poe_all * 100, False_rate_tau_SNRhigh * 100,'r^','Linewidth',lw,'Markersize',msz)
plot(Poe_all * 100, False_rate_CA_SNRhigh * 100,'b^','Linewidth',lw,'Markersize',msz)
legend('${\rm P}_{\rm FA}$ (nominal)', 'NOMP (14 dB)', ...
'NOMP-CFAR (14 dB)', 'NOMP (15 dB)', 'NOMP-CFAR (15 dB)', ...
'NOMP (18 dB)', 'NOMP-CFAR (18 dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('nominal $\bar{\rm P}_{\rm FA}$ (percent)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured $\bar{\rm P}_{\rm FA}$ (percent)', 'Interpreter', 'latex', 'Fontsize', fsz)

figure;
plot(P_FA_all, PDveclow_CFAR, 'k', 'Linewidth', lw);
hold on;
plot(P_FA_all, PDvecmedium_CFAR, 'k--', 'Linewidth', lw);
plot(P_FA_all, PDvechigh_CFAR, 'k:', 'Linewidth', lw);
plot(P_FA_all, Detect_rate_tau_SNRlow, 'ro', 'Linewidth', lw, 'Markersize', msz)
plot(P_FA_all, Detect_rate_CA_SNRlow, 'bo', 'Linewidth', lw, 'Markersize', msz)
plot(P_FA_all, Detect_rate_tau_SNRmedium, 'r+', 'Linewidth', lw, 'Markersize', msz)
plot(P_FA_all, Detect_rate_CA_SNRmedium, 'b+', 'Linewidth', lw, 'Markersize', msz)
plot(P_FA_all, Detect_rate_tau_SNRhigh, 'r^', 'Linewidth', lw, 'Markersize', msz)
plot(P_FA_all, Detect_rate_CA_SNRhigh, 'b^', 'Linewidth', lw, 'Markersize', msz)
legend('Theoretical ($14$ dB)', 'Theoretical ($15$ dB)', ...
'Theoretical ($18$ dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('nominal $\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured $\bar{\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)

% 'NOMP ($14$ dB)', 'NOMP-CFAR ($14$ dB)', ...
% 'NOMP ($15$ dB)', 'NOMP-CFAR ($15$ dB)', ...
% 'NOMP ($18$ dB)', 'NOMP-CFAR ($18$ dB)', ...


% figure;
% plot(P_FA_all, PDvecmedium_tau, '--k', 'Linewidth', lw);
% hold on;
% plot(P_FA_all, Detect_rate_tau_SNRmedium, 'r+', 'Linewidth', lw, 'Markersize', msz)
% plot(P_FA_all, PDvecmedium_CFAR, 'k', 'Linewidth', lw);
% plot(P_FA_all, Detect_rate_CA_SNRmedium, 'b+', 'Linewidth', lw, 'Markersize', msz)
% legend('${\rm P}_{\rm D}$ (NOMP, $15$ dB)', 'NOMP ($15$ dB)', ...
%     '${\rm P}_{\rm D}$ (NOMP-CFAR, $15$ dB)', 'NOMP-CFAR ($15$ dB)', ...
%     'Interpreter', 'latex', 'Fontsize', fsz)
% xlabel('nominal $\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)
% ylabel('measured $\bar{\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)

% figure;
% plot(P_FA_all, PDvechigh_tau, '--k', 'Linewidth', lw);
% hold on;
% plot(P_FA_all, Detect_rate_tau_SNRhigh, 'r^', 'Linewidth', lw, 'Markersize', msz)
% plot(P_FA_all, PDvechigh_CFAR, 'k', 'Linewidth', lw);
% plot(P_FA_all, Detect_rate_CA_SNRhigh, 'b^', 'Linewidth', lw, 'Markersize', msz)
% legend('${\rm P}_{\rm D}$ (NOMP, $18$ dB)', 'NOMP ($18$ dB)', ...
%     '${\rm P}_{\rm D}$ (NOMP-CFAR, $18$ dB)', 'NOMP-CFAR ($18$ dB)', ...
%     'Interpreter', 'latex', 'Fontsize', fsz)
% xlabel('nominal $\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)
% ylabel('measured $\bar{\rm P}_{\rm D}$', 'Interpreter', 'latex', 'Fontsize', fsz)

