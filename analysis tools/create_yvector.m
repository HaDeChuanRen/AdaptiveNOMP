function [y, omega_true, gain_true] = create_yvector(K, S_snap, SNR, sigma_n, Smat_com)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
    [My, Nx] = size(Smat_com);
    omega_true = zeros(K, 1);
    omega_min = 2 * pi / Nx;


    y_noise = sqrt(sigma_n / 2) * (randn(My, S_snap) + 1j*randn(My, S_snap));

    if K > 0
        omega_true(1) = pi * (2 * rand - 1);
        for k = 2 : K
            th = pi * (2 * rand - 1);
            while min(abs((wrapToPi(th - omega_true(1 : k - 1))))) < omega_min
                th = pi * (2*rand - 1);
            end
            omega_true(k) = th;
        end
        omega_true = wrapTo2Pi(omega_true);
        gain_true = bsxfun(@times, sqrt(sigma_n) * (10 .^ (SNR / 20)),...
        exp(1j*2*pi*rand(K, S_snap)));  % K* S_snap
        y_full = exp(1j * (0:(Nx - 1)).' * omega_true.') / sqrt(Nx) * gain_true;
        y = Smat_com * y_full + y_noise;
    else
        gain_true = [];
        omega_true = [];
        y = y_noise;
    end
end

