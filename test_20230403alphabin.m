% last update: 2023/4/3

clc; clear; close all;
addpath('analysis tools\')
addpath('NOMP1D tools\')
MC = 30000;

N_r = 60;
Pfa = 0.01;
Nx = 256;
alpha_set = alpha_PoebyS(Pfa, Nx, N_r);
bvec = zeros(N_r, 1);
nvec = (1 : N_r)';

for n = 2 : N_r
    bvec(n) = N_r * log((n * alpha_set + N_r) ./ (alpha_set + N_r)) + ...
        sum(log(1 : (N_r - n))) - sum(log((n + 1) : (N_r - 1)));
end

lw = 2;
fsz = 12;
msz = 8;
% figure;
% stem(nvec, bvec, '-b', 'Linewidth', lw, 'Markersize', msz);
% xlabel('$n$', 'Interpreter', 'latex', 'Fontsize', fsz)
% ylabel('$b(n)$', 'Interpreter', 'latex', 'Fontsize', fsz)

sigma_n = 1;

% CFAR detector create
num_train = N_r / 2;
num_guard = 4;
cfar_detector = phased.CFARDetector('NumTrainingCells', N_r, 'NumGuardCells',...
    2 * num_guard);
cfar_detector.Method = 'CA';
cfar_detector.ThresholdFactor = 'Auto';
cfar_detector.ProbabilityFalseAlarm = Pfa / Nx;
cfar_detector.ThresholdOutputPort = true;
cfar_detector.NoisePowerOutputPort = true;

% other NOMP-CFAR parameters
Smat_com = eye(Nx);
K_max = 3;

Falmat_alpha = zeros(MC, 1);
Falmat_FFT = zeros(MC, 1);

tic;
for mc_idx = 1 : MC
    hwaitbar = waitbar(mc_idx / MC);
    y_noise = sqrt(sigma_n / 2) * (randn(Nx, 1) + 1j * randn(Nx, 1));
    y_absfft = abs(fft(y_noise));
    prob_ind_ext = repmat(y_absfft .^ 2, [3, 1]);
    peak_grid = cfar_detector(prob_ind_ext, (Nx + 1 : 2 * Nx)');
    Falmat_FFT(mc_idx) = squeeze(sum(peak_grid));

    [omega_alpha, ~, ~] = MNOMP_CFAR_alpha(y_noise, Smat_com, alpha_set, N_r,...
        K_max);
    Falmat_alpha(mc_idx) = length(omega_alpha);

end
time_MC = toc
delete(hwaitbar);

Pfamea_FFT = mean(Falmat_FFT);
Pfamea_alpha = mean(Falmat_alpha);









