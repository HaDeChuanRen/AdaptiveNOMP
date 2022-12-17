% last update: 2022/12/16

clc; clear; close all;

Nx = 256;
S_snap = 1;
sigma_n = 10;
SNR = 24;
Smat_com = eye(Nx);
K = 16;


y_noise = sqrt(sigma_n / 2) .* (randn(Nx, S_snap) + 1j * randn(Nx, S_snap));
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

figure;
qqplot(real(y_noise))

% figure;
% qqplot(real(y_full))





