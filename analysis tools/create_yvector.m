function [y, omega_true, gain_true] = create_yvector(K, T, N, SNR, sigma_n, S)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    omega_true = zeros(K, 1);
    omega_min = 2 * pi / N;
    omega_true(1) = pi * (2 * rand - 1);
    for k = 2 : K
        th = pi * (2 * rand - 1);
        while min(abs((wrapToPi(th - omega_true(1 : k - 1))))) < omega_min
            th = pi * (2*rand - 1);
        end
        omega_true(k) = th;
    end
    omega_true = wrapTo2Pi(omega_true);
    gain_true = bsxfun(@times, sqrt(sigma_n) * (10.^(SNR / 20)), exp(1j*2*pi*rand(K,T)));  % K*T
    % original signal
    y_full = exp(1j * (0:(N-1)).' * omega_true.')/sqrt(N) * gain_true;
    % y_full = zeros(N, T);
    noise = sqrt(sigma_n / 2) * (randn(N, T) + 1j*randn(N, T));
    y_noisy = S * y_full + noise;
    y = y_noisy;
end

