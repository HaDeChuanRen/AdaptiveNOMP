% last update: 2022/12/16

clc; clear; close all;

Nx = 256;
S_snap = 1;
sigma_n = 10;
SNR = 24;
Smat_com = eye(Nx);
K = 16;


y_noise = sqrt(sigma_n / 2) * (randn(Nx, S_snap) + 1j * randn(Nx, S_snap));
omega_true = zeros(K, 1);
omega_min = 2 * pi / Nx;
omega_true(1) = pi * (2 * rand - 1);
for k = 2 : K
    th = pi * (2 * rand - 1);
    while min(abs((wrapToPi(th - omega_true(1 : k - 1))))) < omega_min
        th = pi * (2 * rand - 1);
    end
    omega_true(k) = th;
end
omega_true = wrapTo2Pi(omega_true);
gain_true = bsxfun(@times, sqrt(sigma_n) * (10 .^ (SNR / 20)), ...
    exp(1j*2*pi*rand(K, S_snap)));  % K * S_snap
y_full = exp(1j * (0:(Nx - 1)).' * omega_true.') / sqrt(Nx) * gain_true;
y = Smat_com * y_full + y_noise;

y_fft = fft(y);

figure;
qqplot(real(y_fft))
title('the qqplot of $\tilde{\mathbf{y}}$', 'Interpreter','latex')
xlabel('')
ylabel('')

figure;
qqplot(real(y_noise))
title('the qqplot of noise')
% xlabel('')
% ylabel('')

% figure;
% qqplot(real(y_full))

K_max = 32;
overSamplingRate = 4;
R_c = 1;
R_s = 3;
sampledManifold = preProcessMeasMat(Smat_com, overSamplingRate);
omegaList = [];
gainList  = [];
y_r = y;
residueList = [y_r(:)' * y_r(:)];

figure;
for cnt_k = 1 : K
    [omega_new, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);

    for i = 1:R_s
        [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
            gain_new, Smat_com, sampledManifold.ant_idx, true);
    end

    omegaList = [omegaList; omega_new]; %#ok<*AGROW>
    gainList  = [gainList; gain_new];

    [omegaList, ~, ~] = refineAll(y_r, omegaList, gainList, Smat_com,...
        sampledManifold.ant_idx, R_s, R_c);

    [omegaList, gainList, y_r, A] = solveLeastSquares(y, omegaList, ...
        Smat_com, sampledManifold.ant_idx);

    residue_new = y_r(:)' * y_r(:);
    residueList = [residueList; residue_new];

    yr_fft = fft(y_r);
    title_fig = ['$\hat{K}$ = ', num2str(cnt_k)];
    subplot(sqrt(K), sqrt(K), cnt_k);
    qqplot(real(yr_fft));
    title(title_fig, 'Interpreter', 'latex');
    xlabel('')
    ylabel('')
end

figure;
for cnt_k = 1 : K
    [omega_new, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);

    for i = 1 : R_s
        [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
            gain_new, Smat_com, sampledManifold.ant_idx, true);
    end

    omegaList = [omegaList; omega_new]; %#ok<*AGROW>
    gainList  = [gainList; gain_new];

    [omegaList, ~, ~] = refineAll(y_r, omegaList, gainList, Smat_com,...
        sampledManifold.ant_idx, R_s, R_c);

    [omegaList, gainList, y_r, A] = solveLeastSquares(y, omegaList, ...
        Smat_com, sampledManifold.ant_idx);

    residue_new = y_r(:)' * y_r(:);
    residueList = [residueList; residue_new];

    yr_fft = fft(y_r);
    title_fig = ['$\hat{K}$ = ', num2str(K + cnt_k)];
    subplot(sqrt(K), sqrt(K), cnt_k);
    qqplot(real(yr_fft));
    title(title_fig, 'Interpreter', 'Latex');
    xlabel('')
    ylabel('')
end

