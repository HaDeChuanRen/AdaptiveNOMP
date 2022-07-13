clc; clear; close all;

load("20220712T142752_mc5000_PDvsSNR.mat")
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


lw = 1.6;
fsz = 12;
msz = 10;

figure(1);
plot(SNR_min_all, P_oe001 * ones(1, Num_SNR), '--k','Linewidth',lw)
hold on;
plot(SNR_min_all, Overest_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, Overest_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
legend('${\rm P}_{\rm OE} = 0.01$', 'NOMP', ...
    'NOMP-CA $(N_r = 50)$', 'Interpreter', 'latex', 'Fontsize', fsz)
xlabel('${\rm SNR}$ (dB)', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('measured ${\rm P}_{\rm OE}$', 'Interpreter', 'latex', 'Fontsize', fsz)


figure(2);
plot(SNR_min_all, Equal_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
hold on;
plot(SNR_min_all, Equal_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('measured ${\rm P}(\hat{K} = K)$','Interpreter','latex','Fontsize',fsz)

figure(3);
plot(SNR_min_all, Detect_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
hold on;
plot(SNR_min_all, Detect_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('measured ${\rm P}_{\rm D}$','Interpreter','latex','Fontsize',fsz)

figure(4);
plot(SNR_min_all, reconErr_tau_001Poe_ave,'-ro','Linewidth',lw,'Markersize',msz)
hold on;
plot(SNR_min_all, reconErr_CA25_001Poe_ave,'-b+','Linewidth',lw,'Markersize',msz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('RMSE($\|\hat{\mathbf{y}} - \mathbf{y}\|$) (dB)','Interpreter','latex','Fontsize',fsz)

% figure(4);
% plot(SNR_min_all, False_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
% hold on;
% plot(SNR_min_all, False_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
% xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
% ylabel('measured ${\rm P}_{\rm FA}$','Interpreter','latex','Fontsize',fsz)

figure(5);
plot(SNR_min_all, CRB_omega,'--k','Linewidth',lw)
hold on;
plot(SNR_min_all, MSEdB_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
plot(SNR_min_all, MSEdB_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
ylabel('MSE($\omega$) (dB)','Interpreter','latex','Fontsize',fsz)

% figure(6);
% plot(SNR_min_all, Miss_rate_tau_001Poe,'-ro','Linewidth',lw,'Markersize',msz)
% hold on;
% plot(SNR_min_all, Miss_rate_CA25_001Poe,'-b+','Linewidth',lw,'Markersize',msz)
% xlabel('${\rm SNR}$ (dB)','Interpreter','latex','Fontsize',fsz)
% ylabel('measured ${\rm P}_{\rm M}$','Interpreter','latex','Fontsize',fsz)




