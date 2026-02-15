% last update: 2023/4/3

clc; clear; close all;
addpath('analysis tools\')
addpath('NOMP1D tools\')
rng(1);
MC = 30000;

Nx = 256;
N_r = 60;
Pfa_all = [1e-3, 2e-3, 3e-3, 5e-3, 8e-3, 1e-2, 1.1e-2, 1.2e-2, 1.3e-2, ...
    1.5e-2, 1.7e-2, 2e-2];
length_Pfa = length(Pfa_all);

bmat = zeros(length_Pfa, N_r - 1);
nvec = (1 : N_r)';

for sp_idx = 1 : length_Pfa
    Pfa = Pfa_all(sp_idx);
    alpha_set = alpha_PoebyS(Pfa, Nx, N_r);
    for n = 1 : N_r - 1
        bmat(sp_idx, n) = - (N_r * log(((n + 1) * alpha_set + N_r) ./ ...
            (alpha_set + N_r)) + sum(log(1 : (N_r - n - 1))) - ...
            sum(log((n + 2) : (N_r - 1))));
    end
end

bvec = max(bmat, [], 2);

lw = 2;
fsz = 12;
msz = 8;
figure;
stem(Pfa_all, bvec, '-b', 'Linewidth', lw, 'Markersize', msz);
xlabel('$\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('$b$', 'Interpreter', 'latex', 'Fontsize', fsz)

sigma_n = 1;

alphavec_exa = zeros(length_Pfa, 1);
alphavec_app = zeros(length_Pfa, 1);

for sp_idx = 1 : length_Pfa
    Pfa = Pfa_all(sp_idx);
    alphavec_exa(sp_idx) = alpha_PoebyS(Pfa, Nx, N_r);
    alphavec_app(sp_idx) = N_r * ((Pfa / Nx) ^ (- 1 / N_r) - 1);
end

figure;
hold on;
plot(Pfa_all, alphavec_app, '--k+', 'Linewidth', lw, 'Markersize', msz)
plot(Pfa_all, alphavec_exa, '-bo', 'Linewidth', lw, 'Markersize', msz)
legend('approximation by eq.(33)', 'exact result by eq.(18)', 'Fontsize', fsz)
xlabel('$\bar{\rm P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('$\alpha$', 'Interpreter', 'latex', 'Fontsize', fsz)


% other NOMP-CFAR parameters
% Smat_com = eye(Nx);
% K_max = 2;

% Falmat_alpha = zeros(MC, length_Pfa);
% Falmat_FFT = zeros(MC, length_Pfa);
% num_train = N_r / 2;
% num_guard = 4;

% tic;
% for sp_idx = 1 : length_Pfa
%     % CFAR detector create
%     Pfa = Pfa_all(sp_idx);
%     alpha_set = alpha_PoebyS(Pfa, Nx, N_r);
%     cfar_detector = phased.CFARDetector('NumTrainingCells', N_r, ...
%         'NumGuardCells', 2 * num_guard);
%     cfar_detector.Method = 'CA';
%     cfar_detector.ThresholdFactor = 'Auto';
%     cfar_detector.ProbabilityFalseAlarm = Pfa / Nx;
%     cfar_detector.ThresholdOutputPort = true;
%     cfar_detector.NoisePowerOutputPort = true;

%     for mc_idx = 1 : MC
%         hwaitbar = waitbar((MC * (sp_idx - 1) + mc_idx) / (MC * length_Pfa));
%         y_noise = sqrt(sigma_n / 2) * (randn(Nx, 1) + 1j * randn(Nx, 1));
%         y_absfft = abs(fft(y_noise));
%         prob_ind_ext = repmat(y_absfft .^ 2, [3, 1]);
%         peak_grid = cfar_detector(prob_ind_ext, (Nx + 1 : 2 * Nx)');
%         Falmat_FFT(mc_idx, sp_idx) = squeeze(sum(peak_grid));

%         [omega_alpha, ~, ~] = MNOMP_CFAR_alpha(y_noise, Smat_com, alpha_set, ...
%             N_r, K_max);
%         Falmat_alpha(mc_idx, sp_idx) = length(omega_alpha);

%     end
% end
% time_MC = toc
% delete(hwaitbar);

% if time_MC > 600 && MC > 1000
%     filename_now = [datestr(now, 30), '_mc', num2str(MC), '_PFAbinvsPFA.mat'];
%     save(filename_now, 'Nx', 'Pfa_all', 'length_Pfa', 'MC', 'Falmat_FFT', ...
%         'Falmat_alpha');
% end

% Pfamea_FFT = mean(Falmat_FFT);
% Pfamea_alpha = mean(Falmat_alpha);

% lw = 2;
% fsz = 12;
% msz = 8;

% figure;
% hold on;
% plot(Pfa_all, Pfa_all, '--k', 'Linewidth', lw, 'Markersize', msz)
% plot(Pfa_all, Pfamea_FFT, '-b+', 'Linewidth', lw, 'Markersize', msz);
% plot(Pfa_all, Pfamea_alpha, '-ro', 'Linewidth', lw, 'Markersize', msz);
% legend('$\bar{\rm P}_{\rm FA}$ (nominal)', 'FFT', 'NOMP-CFAR', ...
%     'Interpreter', 'latex', 'Fontsize', fsz)
% xlabel('$\bar{\rm P}_{\rm FA}$ (nominal)', 'Interpreter', 'latex', ...
%     'Fontsize', fsz)
% ylabel('$\bar{\rm P}_{\rm FA}$ (measured)', 'Interpreter', 'latex', ...
%     'Fontsize', fsz)




